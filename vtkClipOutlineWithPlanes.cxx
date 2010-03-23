/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkClipOutlineWithPlanes.cxx,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkClipOutlineWithPlanes.h"

#include "vtkDataSet.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkImageData.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkUnsignedCharArray.h"
#include "vtkDoubleArray.h"
#include "vtkPlaneCollection.h"
#include "vtkMath.h"
#include "vtkMergePoints.h"
#include "vtkGenericCell.h"
#include "vtkPolygon.h"
#include "vtkLine.h"
#include "vtkMatrix4x4.h"
#include "vtkDelaunay2D.h"
//#include "vtkXMLPolyDataWriter.h"

#include "vtkstd/vector"
#include "vtkstd/algorithm"

vtkCxxRevisionMacro(vtkClipOutlineWithPlanes, "$Revision: 1.9 $");
vtkStandardNewMacro(vtkClipOutlineWithPlanes);

vtkCxxSetObjectMacro(vtkClipOutlineWithPlanes,ClippingPlanes,vtkPlaneCollection);

//----------------------------------------------------------------------------
vtkClipOutlineWithPlanes::vtkClipOutlineWithPlanes()
{
  this->ClippingPlanes = 0;
  this->GenerateScalars = 0;
  this->GenerateFaces = 0;
  this->ActivePlaneId = -1;

  this->BaseColor[0] = 1.0;
  this->BaseColor[1] = 0.0;
  this->BaseColor[2] = 0.0;

  this->ClipColor[0] = 1.0;
  this->ClipColor[1] = 0.5;
  this->ClipColor[2] = 0.0;

  this->ActivePlaneColor[0] = 1.0;
  this->ActivePlaneColor[1] = 1.0;
  this->ActivePlaneColor[2] = 0.0;

  // A whole bunch of objects needed during execution
  this->Locator = 0;
  this->CellClipScalars = 0;
  this->IdList = 0;
  this->CellArray = 0;
  this->Polygon = 0;
  this->Cell = 0;
  this->Delaunay = 0;
}

//----------------------------------------------------------------------------
vtkClipOutlineWithPlanes::~vtkClipOutlineWithPlanes()
{
  if (this->ClippingPlanes) { this->ClippingPlanes->Delete(); } 

  if (this->Locator) { this->Locator->Delete(); }
  if (this->CellClipScalars) { this->CellClipScalars->Delete(); }
  if (this->IdList) { this->IdList->Delete(); }
  if (this->CellArray) { this->CellArray->Delete(); }
  if (this->Polygon) { this->Polygon->Delete(); }
  if (this->Cell) { this->Cell->Delete(); }
  if (this->Delaunay) { this->Delaunay->Delete(); }
}

//----------------------------------------------------------------------------
void vtkClipOutlineWithPlanes::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "ClippingPlanes: ";
  if (this->ClippingPlanes)
    {
    os << this->ClippingPlanes << "\n";
    }
  else
    {
    os << "(none)\n";
    } 

  os << indent << "GenerateFaces: "
     << (this->GenerateFaces ? "On\n" : "Off\n" );

  os << indent << "GenerateScalars: "
     << (this->GenerateScalars ? "On\n" : "Off\n" );

  os << indent << "BaseColor: " << this->BaseColor[0] << ", "
     << this->BaseColor[1] << ", " << this->BaseColor[2] << "\n";

  os << indent << "ClipColor: " << this->ClipColor[0] << ", "
     << this->ClipColor[1] << ", " << this->ClipColor[2] << "\n";

  os << indent << "ActivePlaneId: " << this->ActivePlaneId << "\n";

  os << indent << "ActivePlaneColor: " << this->ActivePlaneColor[0] << ", "
     << this->ActivePlaneColor[1] << ", " << this->ActivePlaneColor[2] << "\n";
}

//----------------------------------------------------------------------------
int vtkClipOutlineWithPlanes::ComputePipelineMTime(
  vtkInformation* vtkNotUsed(request),
  vtkInformationVector** vtkNotUsed(inputVector),
  vtkInformationVector* vtkNotUsed(outputVector),
  int vtkNotUsed(requestFromOutputPort),
  unsigned long* mtime)
{
  unsigned long mTime = this->GetMTime();

  vtkPlaneCollection *planes = this->ClippingPlanes;
  vtkPlane *plane = 0;

  if (planes)
    {
    unsigned long planesMTime = planes->GetMTime();
    if (planesMTime > mTime)
      {
      mTime = planesMTime;
      }

    vtkCollectionSimpleIterator iter;
    planes->InitTraversal(iter);
    while ( (plane = planes->GetNextPlane(iter)) )
      {
      unsigned long planeMTime = plane->GetMTime();
      if (planeMTime > mTime)
        {
        mTime = planeMTime;
        }
      }
    }

  *mtime = mTime;

  return 1;
}

//----------------------------------------------------------------------------
int vtkClipOutlineWithPlanes::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // Get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // Get the input and output
  vtkPolyData *input = vtkPolyData::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData *output = vtkPolyData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  // Compute the tolerance based on the data bounds
  double bounds[6];
  input->GetBounds(bounds);
  double tol = 0;
  for (int dim = 0; dim < 3; dim++)
    {
    double d = bounds[2*dim + 1] - bounds[2*dim];
    tol += d*d;
    }
  tol = sqrt(tol)*1e-5;

  // If the input has no line geometry, there's nothing to do
  if (!input->GetLines())
    {
    output->SetPoints(0);
    output->SetLines(0);
    output->SetPolys(0);
    output->GetCellData()->SetScalars(0);

    return 1;
    }

  // The points, forced to double precision.
  vtkPoints *points = vtkPoints::New();
  points->SetDataTypeToDouble();
  vtkPoints *inputPoints = input->GetPoints();
  vtkIdType numPts = inputPoints->GetNumberOfPoints();
  points->SetNumberOfPoints(numPts);
  for (vtkIdType ptId = 0; ptId < numPts; ptId++)
    {
    double point[3];
    inputPoints->GetPoint(ptId, point);
    points->SetPoint(ptId, point);
    } 

  // The cell scalars
  vtkUnsignedCharArray *lineScalars = 0;
  vtkUnsignedCharArray *polyScalars = 0;
  vtkUnsignedCharArray *inputScalars = 0;

  // For input scalars: the offsets
  vtkIdType firstLineScalar = 0;
  vtkIdType firstPolyScalar = 0;

  // Make the colors to be used on the data.
  unsigned char colors[3][3];
  this->CreateColorValues(this->BaseColor, this->ClipColor,
                          this->ActivePlaneColor, colors);

  // This is set if we have to work with scalars.  The input scalars
  // will be utilized if they are unsigned char with 3 components
  if (this->GenerateScalars)
    {
    // Make the scalars
    lineScalars = vtkUnsignedCharArray::New();
    lineScalars->SetNumberOfComponents(3);

    vtkDataArray *tryInputScalars = input->GetCellData()->GetScalars();
    // Get input scalars if they are RGB color scalars
    if (tryInputScalars && tryInputScalars->IsA("vtkUnsignedCharArray") &&
        tryInputScalars->GetNumberOfComponents() == 3)
      {
      inputScalars = static_cast<vtkUnsignedCharArray *>(
        input->GetCellData()->GetScalars());

      vtkIdType numVerts = 0;
      vtkIdType numLines = 0;
      vtkCellArray *tmpCellArray = 0;
      if ( (tmpCellArray = input->GetVerts()) )
        {
        numVerts = tmpCellArray->GetNumberOfCells();
        }
      if ( (tmpCellArray = input->GetLines()) )
        {
        numLines = tmpCellArray->GetNumberOfCells();
        }
      firstLineScalar = numVerts;
      firstPolyScalar = numVerts + numLines;
      }
    }

  // Break the input lines into segments, generate scalars for lines
  vtkCellArray *lines = vtkCellArray::New();
  this->BreakPolylines(input->GetLines(), lines, inputScalars,
                       firstLineScalar, lineScalars, colors[0]);

  // Copy the optional polygons (no triangle strips yet)
  vtkCellArray *polys = 0;
  if (input->GetPolys())
    {
    // If there are line scalars, then poly scalars are needed too
    if (lineScalars)
      {
      polyScalars = vtkUnsignedCharArray::New();
      }

    polys = vtkCellArray::New();
    this->CopyPolygons(input->GetPolys(), polys, inputScalars,
                       firstPolyScalar, polyScalars, colors[0]);
    }

  // Get the clipping planes
  vtkPlaneCollection *planes = this->ClippingPlanes;

  // If there are no planes, then do no clipping
  if (!planes)
    {
    output->SetPoints(points);
    points->Delete();
    output->SetLines(lines);
    lines->Delete();
    output->GetCellData()->SetScalars(lineScalars);
    if (lineScalars) { lineScalars->Delete(); }
    if (polys) { polys->Delete(); }
    if (polyScalars) { polyScalars->Delete(); }

    return 1;
    }

  // The point scalars, for clipping.
  vtkDoubleArray *pointScalars = vtkDoubleArray::New();

  // Arrays for storing the clipped lines and polys.
  vtkCellArray *newLines = vtkCellArray::New();
  vtkCellArray *newPolys = 0;
  if (polys)
    {
    newPolys = vtkCellArray::New();
    }

  // Make the locator and the points
  if (this->Locator == 0) { this->Locator = vtkMergePoints::New(); }
  vtkIncrementalPointLocator *locator = this->Locator;
  vtkPoints *newPoints = vtkPoints::New();
  newPoints->SetDataTypeToDouble();

  // These hold the point and cell scalars, they're needed for clipping
  vtkPointData *inPointData = vtkPointData::New();
  inPointData->CopyScalarsOn();
  inPointData->SetScalars(pointScalars);
  pointScalars->Delete();
  pointScalars = 0;

  vtkCellData *inLineData = vtkCellData::New();
  inLineData->CopyScalarsOn();
  inLineData->SetScalars(lineScalars);
  if (lineScalars)
    {
    lineScalars->Delete();
    lineScalars = 0;
    }

  vtkCellData *inPolyData = vtkCellData::New();
  inPolyData->CopyScalarsOn();
  inPolyData->SetScalars(polyScalars);
  if (polyScalars)
    {
    polyScalars->Delete();
    polyScalars = 0;
    }

  vtkPointData *outPointData = vtkPointData::New();
  outPointData->CopyScalarsOn();

  vtkCellData *outLineData = vtkCellData::New();
  outLineData->CopyScalarsOn();

  vtkCellData *outPolyData = vtkCellData::New();
  outPolyData->CopyScalarsOn();

  // Go through the planes and clip the lines with each plane
  vtkCollectionSimpleIterator iter;
  planes->InitTraversal(iter);
  vtkPlane *plane;
  for (int planeId = 0; (plane = planes->GetNextPlane(iter)); planeId++)
    {
    // Is this the active plane?
    int active = (planeId == this->ActivePlaneId);

    // Convert the plane into an easy-to-evaluate function
    double pc[4];
    plane->GetNormal(pc);
    pc[3] = -vtkMath::Dot(pc, plane->GetOrigin());

    // Create the clip scalars by evaluating the plane at each point
    vtkIdType numPoints = points->GetNumberOfPoints();
    pointScalars = static_cast<vtkDoubleArray *>(inPointData->GetScalars());
    pointScalars->SetNumberOfValues(numPoints);
    for (vtkIdType pointId = 0; pointId < numPoints; pointId++)
      {
      double p[3];
      points->GetPoint(pointId, p);
      double val = p[0]*pc[0] + p[1]*pc[1] + p[2]*pc[2] + pc[3];
      pointScalars->SetValue(pointId, val);
      }

    // Prepare the locator
    locator->InitPointInsertion(newPoints, input->GetBounds());

    // Prepare the output scalars
    outPointData->InterpolateAllocate(inPointData, 0, 0);
    outLineData->CopyAllocate(inLineData, 0, 0);
    outPolyData->CopyAllocate(inPolyData, 0, 0);

    // Clip the lines
    this->ClipAndContourCells(
      points, pointScalars, locator, 1, lines, 0, newLines,
      inPointData, outPointData, inLineData, 0, outLineData);

    if (polys)
      {
      // Get the number of lines remaining after the clipping
      vtkIdType numClipLines = newLines->GetNumberOfCells();

      // Cut the polys to generate more lines
      this->ClipAndContourCells(
        points, pointScalars, locator, 2, polys, newPolys, newLines,
        inPointData, outPointData, inPolyData, outPolyData, outLineData);
      
      // Set new scalars for the contour lines
      vtkUnsignedCharArray *scalars =
        vtkUnsignedCharArray::SafeDownCast(outLineData->GetScalars());

      if (scalars)
        {
        // Set the color to the active color if plane is active
        unsigned char *color = colors[1+active];
        unsigned char *activeColor = colors[2];

        vtkIdType numLines = newLines->GetNumberOfCells();
        for (vtkIdType lineId = numClipLines; lineId < numLines; lineId++)
          {
          unsigned char oldColor[3];
          scalars->GetTupleValue(lineId, oldColor);
          if (oldColor[0] != activeColor[0] ||
              oldColor[1] != activeColor[1] ||
              oldColor[2] != activeColor[2])
            {
            scalars->SetTupleValue(lineId, color);
            }
          }
        }

      // Generate new polys from the cut lines
      this->MakeCutPolys(newPoints, newLines, numClipLines, newPolys, pc,
                         outPolyData, colors[1 + active]);

      }

    // Swap the lines, points, etcetera: old output becomes new input
    vtkPoints *tmp0 = points;
    points = newPoints;
    newPoints = tmp0;
    newPoints->Initialize();

    vtkCellArray *tmp1 = lines;
    lines = newLines;
    newLines = tmp1;
    newLines->Initialize();

    if (polys)
      {
      vtkCellArray *tmp2 = polys;
      polys = newPolys;
      newPolys = tmp2;
      newPolys->Initialize();
      }

    vtkPointData *tmp3 = inPointData;
    inPointData = outPointData;
    outPointData = tmp3;
    outPointData->Initialize();

    vtkCellData *tmp4 = inLineData;
    inLineData = outLineData;
    outLineData = tmp4;
    outLineData->Initialize();

    vtkCellData *tmp5 = inPolyData;
    inPolyData = outPolyData;
    outPolyData = tmp5;
    outPolyData->Initialize();
    }

  output->SetPoints(points);
  points->Delete();

  output->SetLines(lines);
  lines->Delete();

  vtkUnsignedCharArray *scalars = 
    vtkUnsignedCharArray::SafeDownCast(inLineData->GetScalars());

  if (this->GenerateFaces)
    {
    output->SetPolys(polys);

    if (polys && scalars)
      {
      vtkUnsignedCharArray *pScalars = 
        vtkUnsignedCharArray::SafeDownCast(inPolyData->GetScalars());

      vtkIdType n = pScalars->GetNumberOfTuples();
      for (vtkIdType i = 0; i < n; i++)
        {
        unsigned char color[3];
        pScalars->GetTupleValue(i, color);
        scalars->InsertNextTupleValue(color);
        }
      }
   }

  if (polys)
    {
    polys->Delete();
    }

  output->GetCellData()->SetScalars(scalars);

  locator->Initialize();
  newPoints->Delete();
  newLines->Delete();
  if (newPolys)
    {
    newPolys->Delete();
    }

  inPointData->Delete();
  outPointData->Delete();
  inLineData->Delete();
  outLineData->Delete();
  inPolyData->Delete();
  outPolyData->Delete();

  return 1;
}

//----------------------------------------------------------------------------
void vtkClipOutlineWithPlanes::CreateColorValues(
  double color1[3], double color2[3], double color3[3],
  unsigned char colors[3][3])
{
  // Convert colors from "double" to "unsigned char"

  double *dcolors[3];
  dcolors[0] = color1;
  dcolors[1] = color2;
  dcolors[2] = color3;

  for (int i = 0; i < 3; i++)
    {
    for (int j = 0; j < 3; j++)
      {
      double val = dcolors[i][j];
      if (val < 0) { val = 0; }
      if (val > 1) { val = 1; }
      colors[i][j] = static_cast<unsigned char>(val*255);
      }
    }
}

//----------------------------------------------------------------------------
void vtkClipOutlineWithPlanes::ClipAndContourCells(
  vtkPoints *points, vtkDoubleArray *pointScalars,
  vtkIncrementalPointLocator *locator, int dimensionality,
  vtkCellArray *inputCells,
  vtkCellArray *outputPolys, vtkCellArray *outputLines, 
  vtkPointData *inPointData, vtkPointData *outPointData,
  vtkCellData *inCellData,
  vtkCellData *outPolyData, vtkCellData *outLineData)
{
  if (this->CellClipScalars == 0) { this->CellClipScalars = vtkDoubleArray::New(); }
  vtkDoubleArray *cellClipScalars = this->CellClipScalars;

  if (this->Cell == 0) { this->Cell = vtkGenericCell::New(); }
  vtkGenericCell *cell = this->Cell;

  if (this->CellArray == 0) { this->CellArray = vtkCellArray::New(); }
  vtkCellArray *outputVerts = this->CellArray;

  vtkCellData *outCellData = outLineData;
  vtkCellArray *outputCells = outputLines;
  if (dimensionality == 2)
    {
    outCellData = outPolyData;
    outputCells = outputPolys;
    }

  vtkIdType numCells = inputCells->GetNumberOfCells();
  inputCells->InitTraversal();
  for (vtkIdType cellId = 0; cellId < numCells; cellId++)
    {
    vtkIdType numPts, *pts;
    inputCells->GetNextCell(numPts, pts);

    // Set the cell type from the dimensionality
    if (dimensionality == 2)
      {
      if (numPts == 3) { cell->SetCellTypeToTriangle(); }
      else if (numPts == 4) { cell->SetCellTypeToQuad(); }
      else { cell->SetCellTypeToPolygon(); cerr << "polygon\n"; }
      }
    else // (dimensionality == 1)
      {
      if (numPts == 2) { cell->SetCellTypeToLine(); }
      else { cell->SetCellTypeToPolyLine(); }
      }

    vtkPoints *cellPts = cell->GetPoints();
    vtkIdList *cellIds = cell->GetPointIds();

    cellPts->SetNumberOfPoints(numPts);
    cellIds->SetNumberOfIds(numPts);
    cellClipScalars->SetNumberOfValues(numPts);

    // Copy everything over to the temporary cell
    for (vtkIdType i = 0; i < numPts; i++)
      {
      double point[3];
      points->GetPoint(pts[i], point);
      cellPts->SetPoint(i, point);
      cellIds->SetId(i, pts[i]);
      double s = pointScalars->GetValue(cellIds->GetId(i));
      cellClipScalars->SetValue(i, s);
      }

    cell->Clip(0, cellClipScalars, locator, outputCells,
               inPointData, outPointData,
               inCellData, cellId, outCellData, 0);

    if (dimensionality == 2)
      {
      cell->Contour(0, cellClipScalars, locator,
                    outputVerts, outputLines, 0,
                    inPointData, outPointData,
                    inCellData, cellId, outLineData);
      }
    }
}             

//----------------------------------------------------------------------------
void vtkClipOutlineWithPlanes::BreakPolylines(
  vtkCellArray *inputLines, vtkCellArray *lines,
  vtkUnsignedCharArray *inputScalars, vtkIdType firstLineScalar,
  vtkUnsignedCharArray *scalars, unsigned char *color)
{
  // The color for the lines
  unsigned char cellColor[3];
  cellColor[0] = color[0];
  cellColor[1] = color[1];
  cellColor[2] = color[2];

  // Break the input lines into segments
  inputLines->InitTraversal();
  vtkIdType cellId = 0;
  vtkIdType npts, *pts;
  while (inputLines->GetNextCell(npts, pts))
    {
    if (inputScalars)
      {
      inputScalars->GetTupleValue(firstLineScalar + cellId++, cellColor);
      }

    for (vtkIdType i = 1; i < npts; i++)
      {
      lines->InsertNextCell(2);
      lines->InsertCellPoint(pts[i-1]);
      lines->InsertCellPoint(pts[i]);

      if (scalars)
        {
        scalars->InsertNextTupleValue(cellColor);
        }
      }
    }
}

//----------------------------------------------------------------------------
void vtkClipOutlineWithPlanes::CopyPolygons(
  vtkCellArray *inputPolys, vtkCellArray *polys,
  vtkUnsignedCharArray *inputScalars, vtkIdType firstPolyScalar,
  vtkUnsignedCharArray *polyScalars, unsigned char color[3])
{
  polys->DeepCopy(inputPolys);

  // If there are line scalars, we also need poly scalars
  if (polyScalars)
    {
    unsigned char scalarValue[3];
    scalarValue[0] = color[0];
    scalarValue[1] = color[1];
    scalarValue[2] = color[2];

    vtkIdType n = polys->GetNumberOfCells();
    polyScalars->SetNumberOfComponents(3);
    polyScalars->SetNumberOfTuples(n);

    if (inputScalars)
      {
      for (vtkIdType i = 0; i < n; i++)
        {
        inputScalars->GetTupleValue(i + firstPolyScalar, scalarValue);
        polyScalars->SetTupleValue(i, scalarValue);
        }
      }
    else
      {
      for (vtkIdType i = 0; i < n; i++)
        {
        polyScalars->SetTupleValue(i, scalarValue);
        }
      }
    }
}

//----------------------------------------------------------------------------
// A helper class: a bitfield that is always as large as needed.
// For our purposes this is much more convenient than a bool vector,
// which would have to be resized and range-checked externally.

class vtkClipOutlineBitArray
{
public:
  void set(size_t bit, int val) {
    size_t n = (bit >> 5);
    size_t i = (bit & 0x1f);
    if (n >= bitstorage.size()) { bitstorage.resize(n+1); }
    unsigned int chunk = bitstorage[n];
    int bitval = 1;
    bitval <<= i;
    if (val) { chunk = chunk | bitval; }
    else { chunk = chunk & ~bitval; }
    bitstorage[n] = chunk;
  };

  int get(size_t bit) {
    size_t n = (bit >> 5);
    size_t i = (bit & 0x1f);
    if (n >= bitstorage.size()) { return 0; }
    unsigned int chunk = bitstorage[n];
    return ((chunk >> i) & 1);
  };

  void clear() {
    bitstorage.clear();
  };

  void merge(vtkClipOutlineBitArray &b) {
    if (b.bitstorage.size() > bitstorage.size()) {
      bitstorage.resize(b.bitstorage.size()); }
    for (size_t i = 0; i < b.bitstorage.size(); i++) {
      bitstorage[i] |= b.bitstorage[i]; }
  }; 

private:
  vtkstd::vector<unsigned int> bitstorage;
};

//----------------------------------------------------------------------------
// A simple typedefs for stl-based polygons
typedef vtkstd::vector<vtkIdType> vtkClipOutlinePoly;
typedef vtkstd::vector<size_t> vtkClipOutlinePolyGroup;

//----------------------------------------------------------------------------
// See below for more documentation of these methods.

// Take a set of lines, join them to make polygons
static void vtkClipOutlineMakePolysFromLines(
  vtkCellArray *lines, vtkIdType firstLine, vtkIdType numLines,
  vtkstd::vector<vtkClipOutlinePoly> &newPolys);

// Check for polygons that contain multiple loops, and split the loops
static void vtkClipOutlineUntangleSelfIntersection(
  vtkstd::vector<vtkClipOutlinePoly> &newPolys);

// Remove points that are in the middle of the polygon's sides.
static void vtkClipOutlineFindTrueEdges(
  vtkstd::vector<vtkClipOutlinePoly> &newPolys, vtkPoints *points);

// Make sure that the sense of the polygons matches the given normal
static void vtkClipOutlineCorrectPolygonSense(
  vtkstd::vector<vtkClipOutlinePoly> &newPolys, vtkPoints *points,
  const double normal[3]);

// Check for polys within other polys, i.e. find polys that are holes
static void vtkClipOutlineMakeHoleyPolys(
  vtkstd::vector<vtkClipOutlinePoly> &polys, vtkPoints *points,
  vtkstd::vector<vtkClipOutlinePolyGroup> &polyGroups,
  const double normal[3]);

// Add cuts between outer and inner polys, thereby joining the hole
// with the outside perimeter, in order to eliminate the holes
static void vtkClipOutlineCutHoleyPolys(
  vtkstd::vector<vtkClipOutlinePoly> &polys, vtkPoints *points,
  vtkstd::vector<vtkClipOutlinePolyGroup> &polyGroups,
  const double normal[3]);

//----------------------------------------------------------------------------
void vtkClipOutlineWithPlanes::MakeCutPolys(
  vtkPoints *points, vtkCellArray *lines, vtkIdType firstLine,
  vtkCellArray *polys, double normal[3], vtkCellData *outCD,
  unsigned char color[3])
{
  // Find the number of lines generated by cutting polygons
  vtkIdType numLines = lines->GetNumberOfCells();
  vtkIdType numNewLines = numLines - firstLine;

  // If no cut lines were generated, there's nothing to do
  if (firstLine >= numLines)
    {
    return;
    }

  // Join all the new lines into connected groups, i.e. polygons.
  // If we are lucky these will be simple, convex polygons.  But
  // we can't count on that.

  vtkstd::vector<vtkClipOutlinePoly> newPolys;
  vtkClipOutlineMakePolysFromLines(lines, firstLine, numNewLines, newPolys);

  // Some polys might be self-intersecting.  Split these into multiple
  // polys at the intersection points.

  vtkClipOutlineUntangleSelfIntersection(newPolys);

  // Some points might be in the middle of straight line segments.
  // These points can be removed without changing the shape of the
  // polys, and removing them makes triangulation more stable.

  vtkClipOutlineFindTrueEdges(newPolys, points);

  // Check polygon orientation against our clip plane normal, and
  // reverse if necessary.

  vtkClipOutlineCorrectPolygonSense(newPolys, points, normal);

  // Next we have to check for polygons with holes, i.e. polygons that
  // have other polygons inside.  Each polygon is "grouped" with the
  // polygons that make up its holes.  Start with each polygon in its
  // own group.

  size_t numNewPolys = newPolys.size();
  vtkstd::vector<vtkClipOutlinePolyGroup> polyGroups(numNewPolys);
  for (size_t i = 0; i < numNewPolys; i++)
    {
    polyGroups[i].push_back(i);
    }

  // Note: this function will reverse interior polys.  Polys inside the
  // interior polys are split off into their own groups.

  vtkClipOutlineMakeHoleyPolys(newPolys, points, polyGroups, normal);

  // Now that the polys have been grouped, there are two methods
  // that can be used to triangulate them:
  //
  // 1) Turn each group into a vtkPolyData and apply Delaunay2D.
  //
  // 2) Create cuts between exterior and interior polygons so that
  //    a simple EarCut triangulation can be used.
  //
  // Right now only the Delaunay code is complete, but unfortunately
  // vtkDelaunay2D is not stable.


  // Make cuts to create simple polygons out of the holey polys
  vtkClipOutlineCutHoleyPolys(newPolys, points, polyGroups, normal);

  // ------ Triangulation code ------

  // The Delaunay requires creating a new set of points that
  // are transformed to the 2D plane.  Only the points that are
  // vertices of our polygons can be included, plus any points
  // that are strictly interior to those polygons.
  // The clipping plane normal and position is used to
  // compute the transform.  A mapping from the new pointIds
  // to the original pointIds will be needed in order to
  // merge the triangulation with our output cell array.

  // Instantiate input polydata for triangulation
  vtkPolyData *sourceData = vtkPolyData::New();
  vtkCellArray *sourcePolys = vtkCellArray::New();
  vtkCellArray *outputPolys = 0;
  vtkPoints *sourcePoints = vtkPoints::New();
  sourcePoints->SetDataTypeToDouble();

  // Build a matrix to transform points into the xy plane.  The normal
  // must be the same as the one we used to set the sense of the polys.
  double matrix[16];
  vtkMatrix4x4::Identity(matrix);
  matrix[8] = -normal[0];
  matrix[9] = -normal[1];
  matrix[10] = -normal[2];
  vtkMath::Perpendiculars(&matrix[8], &matrix[0], &matrix[4], 0);

  // The polygons have been grouped prior to triangulation.  A group
  // of simple polygons is used to represent polygons with holes.

  for (size_t groupId = 0; groupId < polyGroups.size(); groupId++)
    {
    if (polyGroups[groupId].size() == 0)
      {
      continue;
      }

    // Make a map for recovering pointIds after triangulation has finished
    vtkstd::vector<vtkIdType> pointMap;
    vtkIdType newPointId = 0;

    sourcePolys->Initialize();
    sourcePoints->Initialize();
    sourceData->Initialize();
    sourceData->SetPoints(sourcePoints);
    sourceData->SetPolys(sourcePolys);

    // Go through all polys, check whether they are in the group.
    for (size_t i = 0; i < polyGroups[groupId].size(); i++)
      {
      size_t polyId = polyGroups[groupId][i];
      size_t n = newPolys[polyId].size();

      // If the poly is a line, then skip it
      if (n < 3) { continue; }

      // Create the source polys for a constrained Delaunay
      sourcePolys->InsertNextCell(n);
      for (size_t j = 0; j < n; j++)
        {
        // Use homogeneous point with Matrix4x4
        double point[4];
        point[3] = 1.0;

        vtkIdType pointId = newPolys[polyId][j];
        points->GetPoint(pointId, point);
        vtkMatrix4x4::MultiplyPoint(matrix, point, point);

        sourcePoints->InsertNextPoint(point);
        sourcePolys->InsertCellPoint(newPointId++);

        // Save the old pointId
        pointMap.push_back(pointId);
        }
      }

    // Just in case there were no good cells
    if (sourcePolys->GetNumberOfCells() == 0)
      {
      continue;
      }

    // If polys has one cell, don't use Delaunay
    if (sourcePolys->GetNumberOfCells() == 1)
      {
      if (sourcePolys->GetNumberOfConnectivityEntries() <= 5)
        {
        // If the cell is a triangle or quad, don't triangulate
        outputPolys = sourcePolys;
        }
      else
        {
        // If it is a single polygon, use ear cut triangulation

        if (this->Polygon == 0) { this->Polygon = vtkPolygon::New(); }
        vtkPolygon *polygon = this->Polygon;

        if (this->IdList == 0) { this->IdList = vtkIdList::New(); }
        vtkIdList *triangles = this->IdList;

        if (this->CellArray == 0) { this->CellArray = vtkCellArray::New(); }
        vtkCellArray *trianglePolys = this->CellArray;

        vtkIdType npts, *pts;
        sourcePolys->GetCell(0, npts, pts);

        polygon->Points->SetDataTypeToDouble();
        polygon->Points->SetNumberOfPoints(npts);
        polygon->PointIds->SetNumberOfIds(npts);

        for (vtkIdType j = 0; j < npts; j++)
          {
          double point[3];
          sourcePoints->GetPoint(pts[j], point);
          polygon->Points->SetPoint(j, point);
          polygon->PointIds->SetId(j, pts[j]);
          }

        triangles->Initialize();
        polygon->Triangulate(triangles);
        vtkIdType m = triangles->GetNumberOfIds();
        trianglePolys->Initialize();
        for (vtkIdType k = 0; k < m; k += 3)
          {
          trianglePolys->InsertNextCell(3);
          trianglePolys->InsertCellPoint(pts[triangles->GetId(k + 0)]);
          trianglePolys->InsertCellPoint(pts[triangles->GetId(k + 1)]);
          trianglePolys->InsertCellPoint(pts[triangles->GetId(k + 2)]);
          }
        outputPolys = trianglePolys;
        }
      }
    else
      {
      cerr << "doing delaunay!!!\n";
      // Delaunay triangulation

      // Debug code: write data sets that cause Delaunay to die
      //static unsigned int counter = 0;
      //cerr << "delaunay! " << counter++ << "\n";
      //vtkXMLPolyDataWriter *writer = vtkXMLPolyDataWriter::New();
      //writer->SetFileName("baddelaunay.vtp");
      //writer->SetInput(sourceData);
      //writer->SetCompressorTypeToNone();
      //writer->SetDataModeToAscii();
      //writer->Modified();
      //writer->Write();
      //writer->Delete();

      if (this->Delaunay == 0) { this->Delaunay = vtkDelaunay2D::New(); }
      vtkDelaunay2D *delaunay = this->Delaunay;
      delaunay->SetTolerance(1e-5);
      delaunay->SetOffset(1.0);
      delaunay->SetInput(sourceData);
      delaunay->SetSource(sourceData);
      delaunay->SetProjectionPlaneMode(VTK_DELAUNAY_XY_PLANE);

      delaunay->Modified();
      delaunay->Update();
      outputPolys = delaunay->GetOutput()->GetPolys();
      }

    // Check to make sure that Delaunay produced something
    if (outputPolys)
      {
      vtkUnsignedCharArray *scalars =
        vtkUnsignedCharArray::SafeDownCast(outCD->GetScalars());

      outputPolys->InitTraversal();
      vtkIdType npts, *pts;
      while (outputPolys->GetNextCell(npts, pts))
        {
        vtkIdType cellId = polys->InsertNextCell(npts);

        for (vtkIdType k = 0; k < npts; k++)
          {
          polys->InsertCellPoint(pointMap[pts[k]]);
          }

        if (scalars)
          {
          scalars->InsertTupleValue(cellId, color);
          }
        }
      }
    }

  // Free up all the memory that we can

  if (this->Polygon)
    {
    this->Polygon->Points->Initialize();
    this->Polygon->PointIds->Initialize();
    }
  if (this->IdList) { this->IdList->Initialize(); }
  if (this->CellArray) { this->CellArray->Initialize(); }

  if (this->Delaunay)
    {
    this->Delaunay->SetInput(0);
    this->Delaunay->SetSource(0);
    this->Delaunay->SetOutput(0);
    }

  sourcePoints->Delete();
  sourcePolys->Delete();
  sourceData->Delete();
}

// ---------------------------------------------------
// Here is the code for creating polygons from line segments.

void vtkClipOutlineMakePolysFromLines(
  vtkCellArray *lines, vtkIdType firstLine, vtkIdType numLines,
  vtkstd::vector<vtkClipOutlinePoly> &newPolys)
{
  // Skip through the cell array until we get to the first line
  lines->InitTraversal();
  vtkIdType npts, *pts;
  for (vtkIdType cellId = 0; cellId < firstLine; cellId++)
    {
    lines->GetNextCell(npts, pts);
    }

  vtkIdType firstLineLoc = lines->GetTraversalLocation();

  // Bitfield for marking lines as used
  vtkClipOutlineBitArray usedLines;

  size_t numNewPolys = 0;
  vtkIdType remainingLines = numLines;
  while (remainingLines > 0)
    {
    // Create a new poly
    newPolys.resize(++numNewPolys);
    vtkClipOutlinePoly &poly = newPolys[numNewPolys-1];

    int completePoly = 0;
    int noLinesMatch = 0;
    while (!completePoly && !noLinesMatch && remainingLines > 0)
      {
      noLinesMatch = 1;
      lines->SetTraversalLocation(firstLineLoc);
      for (vtkIdType lineId = 0; lineId < numLines; lineId++)
        {
        lines->GetNextCell(npts, pts);

        if (usedLines.get(lineId))
          {
          continue;
          }

        // Number of points in the poly
        size_t npoly = poly.size();

        // Other useful counters
        vtkIdType n = npts;
        vtkIdType m = npoly/2;

        int usedLine = 1;

        if (poly.size() == 0)
          {
          for (vtkIdType i = 0; i < npts; i++)
            { 
            poly.push_back(pts[i]);
            }
          }
        else if (pts[0] == poly[npoly-1])
          {
          if (pts[npts-1] == poly[0])
            {
            n = n-1;
            completePoly = 1;
            }
          for (vtkIdType k = 1; k < n; k++)
            {
            poly.push_back(pts[k]);
            }
          }
        else if (pts[npts-1] == poly[npoly-1])
          {
          if (pts[0] == poly[0])
            {
            n = n-1;
            completePoly = 1;
            }
          for (vtkIdType k = n-1; k >= 1; k--)
            {
            poly.push_back(pts[k]);
            }
          }
        else if (pts[0] == poly[0])
          {
          for (vtkIdType j = 0; j < m; j++)
            {
            vtkIdType tmp = poly[j];
            poly[j] = poly[npoly - j - 1];
            poly[npoly - j - 1] = tmp;
            }
          for (vtkIdType k = 1; k < n; k++)
            {
            poly.push_back(pts[k]);
            }
          }
        else if (pts[0] == poly[npoly-1])
          {
          for (vtkIdType j = 0; j < m; j++)
            {
            vtkIdType tmp = poly[j];
            poly[j] = poly[npoly - j - 1];
            poly[npoly - j - 1] = tmp;
            }
          for (vtkIdType k = n-1; k >= 1; k--)
            {
            poly.push_back(pts[k]);
            }
          }
        else
          {
          usedLine = 0;
          }

        if (usedLine)
          {
          noLinesMatch = 0;
          usedLines.set(lineId, 1);
          remainingLines--;
          }
        }
      }
    }
}

// ---------------------------------------------------
// Check for self-intersection. Split the figure-eights.
// This assumes that all intersections occur at existing
// vertices, i.e. no new vertices will be created.

void vtkClipOutlineUntangleSelfIntersection(
  vtkstd::vector<vtkClipOutlinePoly> &newPolys)
{
  size_t numNewPolys = newPolys.size();
  for (size_t i = 0; i < numNewPolys; i++)
    {
    size_t n = newPolys[i].size();

    int foundMatch = 0;
    size_t idx1 = 0;
    size_t idx2 = 0;

    for (idx1 = 0; idx1 < n; idx1++)
      {
      vtkIdType firstId = newPolys[i][idx1];

      for (idx2 = idx1+1; idx2 < n ; idx2++)
        {
        vtkIdType secondId = newPolys[i][idx2];

        if (firstId == secondId)
          {
          foundMatch = 1;
          break;
          }
        }

      if (foundMatch) { break; }
      }

    if (foundMatch)
      {
      // Split off a new poly
      size_t m = idx2 - idx1;

      newPolys.resize(++numNewPolys);
      newPolys[numNewPolys-1].resize(n - m);

      for (size_t j = 0; j < idx1; j++)
        {
        newPolys[numNewPolys-1][j] = newPolys[i][j];
        }
      for (size_t k = idx2; k < n; k++)
        {
        newPolys[numNewPolys-1][k - m] = newPolys[i][k];
        }

      // The current poly, which is now intersection-free
      for (size_t l = 0; l < m; l++)
        {
        newPolys[i][l] = newPolys[i][l + idx1];
        }
      newPolys[i].resize(m);
      }
    }
}

// ---------------------------------------------------
// The polygons might have a lot of extra points, i.e. points
// in the middle of the edges.  Remove those points.

void vtkClipOutlineFindTrueEdges(
  vtkstd::vector<vtkClipOutlinePoly> &newPolys, vtkPoints *points)
{
  // Tolerance^2 for cross product to see if vectors are parallel
  const double tol2 = 1e-10;

  size_t numNewPolys = newPolys.size();
  for (size_t i = 0; i < numNewPolys; i++)
    {
    vtkClipOutlinePoly newPoly;

    size_t n = newPolys[i].size();

    if (n < 3) { continue; }

    double p0[3], p1[3], p2[3];
    double v1[3], v2[3], v[3];
    double l1, l2;

    points->GetPoint(newPolys[i][n-1], p0);
    points->GetPoint(newPolys[i][0], p1);
    v1[0] = p1[0] - p0[0];  v1[1] = p1[1] - p0[1];  v1[2] = p1[2] - p0[2];  
    l1 = vtkMath::Dot(v1, v1);

    for (size_t j = 0; j < n; j++)
      {
      size_t k = j+1;
      if (k == n) { k = 0; }

      points->GetPoint(newPolys[i][k], p2);
      v2[0] = p2[0] - p1[0];  v2[1] = p2[1] - p1[1];  v2[2] = p2[2] - p1[2];  
      l2 = vtkMath::Dot(v2, v2);

      // magnitude of cross product is |v1||v2|sin(angle)
      vtkMath::Cross(v1, v2, v);
      // s2 is square of magnitude of cross product
      double s2 = vtkMath::Dot(v, v);
      // c is |v1||v2|cos(angle)
      double c = vtkMath::Dot(v1, v2);

      // Keep the point if angle is greater than tolerance
      if (c < 0 || s2 > l1*l2*tol2)
        {
        newPoly.push_back(newPolys[i][j]);
        }

      p1[0] = p2[0]; p1[1] = p2[1]; p1[2] = p2[2];
      v1[0] = v2[0]; v1[1] = v2[1]; v1[2] = v2[2];
      l1 = l2;
      }

    newPolys[i] = newPoly;
    }
}

// ---------------------------------------------------
// Correct the sense of the polygons, by making sure that their
// normal matches the given normal.

void vtkClipOutlineCorrectPolygonSense(
  vtkstd::vector<vtkClipOutlinePoly> &newPolys, vtkPoints *points,
  const double normal[3])
{
  size_t numNewPolys = newPolys.size();
  for (size_t i = 0; i < numNewPolys; i++)
    {
    size_t n = newPolys[i].size();

    if (n < 3) { continue; }

    // Compute the normal, reverse polygon if necessary

    double pnormal[3], p0[3], p1[3], p2[3], v1[3], v2[3], v[3];
    pnormal[0] = 0.0; pnormal[1] = 0.0; pnormal[2] = 0.0;

    points->GetPoint(newPolys[i][0], p0);
    points->GetPoint(newPolys[i][1], p1);
    v1[0] = p1[0] - p0[0];  v1[1] = p1[1] - p0[1];  v1[2] = p1[2] - p0[2];  

    for (size_t jj = 2; jj < n; jj++)
      {
      points->GetPoint(newPolys[i][jj], p2);
      v2[0] = p2[0] - p0[0];  v2[1] = p2[1] - p0[1];  v2[2] = p2[2] - p0[2];  
      vtkMath::Cross(v1, v2, v);
      pnormal[0] += v[0]; pnormal[1] += v[1]; pnormal[2] += v[2];
      p1[0] = p2[0]; p1[1] = p2[1]; p1[2] = p2[2];
      v1[0] = v2[0]; v1[1] = v2[1]; v1[2] = v2[2];
      }

    // The cut normal is inward, the poly normal should be outward
    if (vtkMath::Dot(normal, pnormal) > 0)
      {
      // Reverse the polygon
      size_t m = n/2;
      for (size_t kk = 0; kk < m; kk++)
        {
        vtkIdType tmpId = newPolys[i][kk];
        newPolys[i][kk] = newPolys[i][n - kk - 1];
        newPolys[i][n - kk - 1] = tmpId;
        }
      }
    }
}

// ---------------------------------------------------
// Check whether innerPoly is inside outerPoly.
// The normal is needed to verify the polygon orientation.
// The values of pp, bounds, and tol2 must be precomputed
// by calling vtkClipOutlinePrepareForPolyInPoly() on outerPoly.

int vtkClipOutlinePolyInPoly(
  const vtkClipOutlinePoly &outerPoly,
  const vtkClipOutlinePoly &innerPoly,
  vtkPoints *points, const double normal[3],
  const double *pp, const double bounds[6],
  double tol2)
{
  // Find a vertex of poly "j" that isn't on the edge of poly "i".
  // This is necessary or the PointInPolygon might return "true"
  // based only on roundoff error.

  double p[3];
  int allPointsOnEdges = 1;
  size_t n = outerPoly.size();
  size_t m = innerPoly.size();

  for (size_t jj = 0; jj < m; jj++)
    {          
    points->GetPoint(innerPoly[jj], p);

    int pointOnEdge = 0;
    double q1[3], q2[3];
    points->GetPoint(outerPoly[n-1], q1);
    for (size_t ii = 0; ii < n; ii++)
      {
      points->GetPoint(outerPoly[ii], q2);
      double t, dummy[3];
      // This method returns distance squared
      if (vtkLine::DistanceToLine(p, q1, q2, t, dummy) < tol2)
        {
        pointOnEdge = 1;
        break;
        }
      q1[0] = q2[0]; q1[1] = q2[1]; q1[2] = q2[2];
      }
    if (!pointOnEdge)
      {
      allPointsOnEdges = 0;
      break;
      }
    }

  if (allPointsOnEdges)
    {
    return 1;
    }

  // There should also be a check to see if all the verts match.
  // If they do, both polys should be removed.
   
  return vtkPolygon::PointInPolygon(p, n, const_cast<double *>(pp),
    const_cast<double *>(bounds), const_cast<double *>(normal));
}

// ---------------------------------------------------
// Precompute values needed for the PolyInPoly check.
// The values that are returned are as follows:
// pp: an array of the polygon vertices
// bounds: the polygon bounds
// tol2: a tolerance value based on the size of the polygon
// (note: pp must be pre-allocated to the 3*outerPoly.size())

void vtkClipOutlinePrepareForPolyInPoly(
  const vtkClipOutlinePoly &outerPoly, vtkPoints *points,
  double *pp, double bounds[6], double &tol2)
{
  // Use pp to store the polygon vertices
  size_t n = outerPoly.size();
  bounds[0] = bounds[2] = bounds[4] = 0;
  bounds[1] = bounds[3] = bounds[5] = 0;

  // Find the bounding box for the polygon
  for (size_t k = 0; k < n; k++)
    {
    double *p = &pp[3*k];
    points->GetPoint(outerPoly[k], p);
    if (k == 0)
      {
      bounds[0] = p[0]; bounds[1] = p[0];
      bounds[2] = p[1]; bounds[3] = p[1];
      bounds[4] = p[2]; bounds[5] = p[2];
      }
    else
      {
      if (p[0] < bounds[0]) { bounds[0] = p[0]; }
      if (p[0] > bounds[1]) { bounds[1] = p[0]; }
      if (p[1] < bounds[2]) { bounds[2] = p[1]; }
      if (p[1] > bounds[3]) { bounds[3] = p[1]; }
      if (p[2] < bounds[4]) { bounds[4] = p[2]; }
      if (p[2] > bounds[5]) { bounds[5] = p[2]; }
      }
    }

  // Compute a tolerance based on the poly size
  double ps[3];
  ps[0] = bounds[1] - bounds[0];
  ps[1] = bounds[3] - bounds[2];
  ps[2] = bounds[5] - bounds[4];

  // Tolerance is for squared distance
  tol2 = (ps[0]*ps[0] + ps[1]*ps[1] + ps[2]*ps[2])*(1e-5 * 1e-5);
}

// ---------------------------------------------------
// Check for polygons within polygons.  Group the polygons
// if they are within each other.  Reverse the sense of 
// the interior "hole" polygons.  A hole within a hole
// will be reversed twice and will become its own group.

void vtkClipOutlineMakeHoleyPolys(
  vtkstd::vector<vtkClipOutlinePoly> &newPolys, vtkPoints *points,
  vtkstd::vector<vtkClipOutlinePolyGroup> &polyGroups,
  const double normal[3])
{
  size_t numNewPolys = newPolys.size();
  if (numNewPolys <= 1)
    {
    return;
    }

  // Use bit arrays to keep track of inner polys
  vtkClipOutlineBitArray polyReversed;
  vtkClipOutlineBitArray innerPolys;

  // Find the maximum poly size
  size_t nmax = 1;
  for (size_t kk = 0; kk < numNewPolys; kk++)
    {
    size_t n = newPolys[kk].size();
    if (n > nmax) { nmax = n; }
    }

  // These are some values needed for poly-in-poly checks
  double *pp = new double[3*nmax];
  double bounds[6];
  double tol2;

  // Go through all polys
  for (size_t i = 0; i < numNewPolys; i++)
    {
    size_t n = newPolys[i].size();

    if (n < 3) { continue; }

    // Precompute some values needed for poly-in-poly checks
    vtkClipOutlinePrepareForPolyInPoly(newPolys[i], points, pp, bounds, tol2);

    // Look for polygons inside of this one
    for (size_t j = 0; j < numNewPolys; j++)
      {
      size_t m = newPolys[j].size();
      if (j == i || m < 3) { continue; }

      // Make sure polygon i is not in polygon j
      int isInteriorPoly = 0;
      for (size_t k = 1; k < polyGroups[j].size(); k++)
        {
        if (polyGroups[j][k] == i)
          {
          isInteriorPoly = 1;
          break;
          }
        }

      if (isInteriorPoly)
        {
        continue;
        }

      if (vtkClipOutlinePolyInPoly(newPolys[i], newPolys[j], points,
                                  normal, pp, bounds, tol2))
        {
        // Mark the inner poly as reversed
        polyReversed.set(j, !polyReversed.get(j));

        // Add to group
        polyGroups[i].push_back(j);
        }
      }
    }

  delete [] pp;

  for (size_t j = 0; j < numNewPolys; j++)
    {
    // Reverse the interior polys, and remove their groups
    if (polyReversed.get(j))
      {
      size_t m = newPolys[j].size();
      size_t m2 = m/2;
      for (size_t k = 0; k < m2; k++)
        {
        vtkIdType tmpId = newPolys[j][k];
        newPolys[j][k] = newPolys[j][m - k - 1];
        newPolys[j][m - k - 1] = tmpId;
        }

      polyGroups[j].clear();
      }
    // Polys inside the interior polys have their own groups, so remove
    // them from this group
    else if (polyGroups[j].size() > 1)
      {
      // Convert the group into a bit array, to make manipulation easier
      innerPolys.clear();
      for (size_t k = 1; k < polyGroups[j].size(); k++)
        {
        innerPolys.set(polyGroups[j][k], 1);
        }

      // Look for non-reversed polys inside this one
      for (size_t kk = 1; kk < polyGroups[j].size(); kk++)
        {
        // jj is the index of the inner poly
        size_t jj = polyGroups[j][kk];
        // If inner poly is not reversed then
        if (!polyReversed.get(jj))
          {
          // Remove that poly and all polys inside of it from the group
          for (size_t ii = 0; ii < polyGroups[jj].size(); ii++)
            {
            innerPolys.set(polyGroups[jj][ii], 0);
            }
          }
        }

      // Use the bit array to recreate the polyGroup
      polyGroups[j].clear();
      polyGroups[j].push_back(j);
      for (size_t jj = 0; jj < numNewPolys; jj++)
        {
        if (innerPolys.get(jj) != 0)
          {
          polyGroups[j].push_back(jj);
          }
        }
      }
    }
}

// ---------------------------------------------------
// Check line segment with point Ids (i, j) to make sure that it
// doesn't cut through the edges of any polys in the group.
// Return value of zero means check failed and the cut is not
// usable.

int vtkClipOutlineCheckCut(
  vtkstd::vector<vtkClipOutlinePoly> &polys, vtkPoints *points,
  vtkClipOutlinePolyGroup &polyGroup, vtkIdType ptId1, vtkIdType ptId2)
{
  double p1[3], p2[3];

  points->GetPoint(ptId1, p1);
  points->GetPoint(ptId2, p2);

  for (size_t i = 0; i < polyGroup.size(); i++)
    {
    vtkClipOutlinePoly &poly = polys[polyGroup[i]];
    size_t n = poly.size();

    double q1[3], q2[3];
    vtkIdType qtId1 = poly[n-1];
    points->GetPoint(qtId1, q1);

    for (size_t j = 0; j < n; j++)
      {
      vtkIdType qtId2 = poly[j];
      points->GetPoint(qtId2, q2);

      // If lines share an endpoint, they can't intersect,
      // so don't bother with the check.
      if (ptId1 != qtId1 && ptId1 != qtId2 &&
          ptId2 != qtId1 && ptId2 != qtId2)
        {
        double u, v;
        if (vtkLine::Intersection(p1, p2, q1, q2, u, v))
          {
          return 0;
          }
        }

      qtId1 = qtId2;
      q1[0] = q2[0]; q1[1] = q2[1]; q1[2] = q2[2];
      }
    }

  return 1;
}

// ---------------------------------------------------
// Check the quality of a cut between an outer and inner polygon.
// Larger values mean that the cut will produce triangles with
// sharp angles.  The range of values is [-1, 1], where the smallest
// values indicate the highest quality.

double vtkClipOutlineCutQuality(
  vtkClipOutlinePoly &outerPoly, vtkClipOutlinePoly &innerPoly,
  size_t i, size_t j, vtkPoints *points)
{
  size_t n = outerPoly.size();
  size_t m = innerPoly.size();

  size_t a = ((i > 0) ? i-1 : n-1);
  size_t b = ((i < n-1) ? i+1 : 0);

  size_t c = ((j > 0) ? j-1 : m-1);
  size_t d = ((j < m-1) ? j+1 : 0);

  double p0[3], p1[3], p2[3];
  points->GetPoint(outerPoly[i], p1);
  points->GetPoint(innerPoly[j], p2);

  double v1[3], v2[3];
  v1[0] = p2[0] - p1[0]; v1[1] = p2[1] - p1[1]; v1[2] = p2[2] - p1[2];

  double l1 = sqrt(vtkMath::Dot(v1, v1));
  double l2;
  double qmax = -l1;
  double q;

  points->GetPoint(outerPoly[a], p0);
  v2[0] = p0[0] - p1[0]; v2[1] = p0[1] - p1[1]; v2[2] = p0[2] - p1[2];
  l2 = sqrt(vtkMath::Dot(v2, v2));
  if (l2 > 0)
    {
    q = vtkMath::Dot(v1, v2)/l2;
    if (q > qmax) { qmax = q; }
    }

  points->GetPoint(outerPoly[b], p0);
  v2[0] = p0[0] - p1[0]; v2[1] = p0[1] - p1[1]; v2[2] = p0[2] - p1[2];
  l2 = sqrt(vtkMath::Dot(v2, v2));
  if (l2 > 0)
    {
    q = vtkMath::Dot(v1, v2)/l2;
    if (q > qmax) { qmax = q; }
    }

  points->GetPoint(innerPoly[c], p0);
  v2[0] = p0[0] - p2[0]; v2[1] = p0[1] - p2[1]; v2[2] = p0[2] - p2[2];
  l2 = sqrt(vtkMath::Dot(v2, v2));
  if (l2 > 0)
    {
    q = vtkMath::Dot(v1, v2)/l2;
    if (q > qmax) { qmax = q; }
    }

  points->GetPoint(innerPoly[d], p0);
  v2[0] = p0[0] - p2[0]; v2[1] = p0[1] - p2[1]; v2[2] = p0[2] - p2[2];
  l2 = sqrt(vtkMath::Dot(v2, v2));
  if (l2 > 0)
    {
    q = vtkMath::Dot(v1, v2)/l2;
    if (q > qmax) { qmax = q; }
    }

  if (l1 > 0)
    {
    return qmax/l1;
    }

  return 1.0;
}

// ---------------------------------------------------
// A simple struct to hold the endpoints of a cut and
// a quality value for the cut
class vtkClipOutlinePolyCut
{
public:
  vtkClipOutlinePolyCut(double q, size_t p1, size_t p2) :
    Quality(q), OuterIdx(p1), InnerIdx(p2) {};

  int operator<(const vtkClipOutlinePolyCut& c) const {
    return (this->Quality < c.Quality); };

  int operator>(const vtkClipOutlinePolyCut& c) const {
    return (this->Quality > c.Quality); };

  double Quality;
  size_t OuterIdx;
  size_t InnerIdx;
};

// ---------------------------------------------------
// After the holes have been identified, make cuts between the
// outer poly and each hole.  Make two cuts per hole.  The only
// strict requirement is that the cut must not intersect any
// edges, but it's best to make sure that no really sharp angles
// are created.

void vtkClipOutlineCutHoleyPolys(
  vtkstd::vector<vtkClipOutlinePoly> &polys, vtkPoints *points,
  vtkstd::vector<vtkClipOutlinePolyGroup> &polyGroups,
  const double normal[3])
{
  // Go through all groups and cut out the first inner poly that is
  // found.  Every time an inner poly is cut out, the groupId counter
  // is reset because a cutting a poly creates a new group.
  size_t groupId = 0;
  while (groupId < polyGroups.size())
    {
    vtkClipOutlinePolyGroup &polyGroup = polyGroups[groupId];

    // Only need to make a cut if the group size is greater than 1
    if (polyGroup.size() > 1)
      {
      // The first member of the group is the outer poly
      size_t outerPolyId = polyGroup[0];
      vtkClipOutlinePoly &outerPoly = polys[outerPolyId];

      // The second member of the group is the first inner poly
      size_t innerPolyId = polyGroup[1];
      vtkClipOutlinePoly &innerPoly = polys[innerPolyId];

      // Make a container for "cut" candidates
      vtkstd::vector<vtkClipOutlinePolyCut> cuts;

      // Brute-force search for potential cuts
      for (size_t j = 0; j < outerPoly.size(); j++)
        {
        for (size_t k = 0; k < innerPoly.size(); k++)
          {
          if (vtkClipOutlineCheckCut(polys, points, polyGroup,
                                     outerPoly[j], innerPoly[k]))
            {
            double q = vtkClipOutlineCutQuality(outerPoly, innerPoly, j, k,
                                                points);
            cuts.push_back(vtkClipOutlinePolyCut(q, j, k));
            }
          }
        }

      // There must be at least 2 valid cuts
      assert(cuts.size() >= 2);

      // Sort the cuts to find the best one
      vtkstd::sort(cuts.begin(), cuts.end());

      vtkIdType ptId1 = outerPoly[cuts[0].OuterIdx];
      vtkIdType ptId2 = innerPoly[cuts[0].InnerIdx];
      double p1[3], p2[3];
      points->GetPoint(ptId1, p1);
      points->GetPoint(ptId2, p2);

      // Find the second-best cut that doesn't intersect the first
      size_t k;
      for (k = 1; k < cuts.size(); k++)
        {
        vtkIdType qtId1 = outerPoly[cuts[k].OuterIdx];
        vtkIdType qtId2 = innerPoly[cuts[k].InnerIdx];

        // Second cut shouldn't share either endpoint with first cut
        if (ptId1 == qtId1 || ptId2 == qtId2) { continue; }

        double q1[3], q2[3];
        points->GetPoint(qtId1, q1);
        points->GetPoint(qtId2, q2);

        double u, v;
        if (!vtkLine::Intersection(p1, p2, q1, q2, u, v))
          {
          // Found a good second cut
          break;
          }
        }

      // Ensure that a second cut was found
      assert(k < cuts.size());

      // Generate new polys from the cuts
      size_t a = cuts[0].OuterIdx;
      size_t b = cuts[k].OuterIdx;
      size_t c = cuts[k].InnerIdx;
      size_t d = cuts[0].InnerIdx;

      size_t n = outerPoly.size();
      size_t m = innerPoly.size();
      size_t idx;

      // Generate poly1
      vtkClipOutlinePoly poly1;
      for (idx = a;;)
        {
        poly1.push_back(outerPoly[idx]);
        if (idx == b) { break; }
        if (++idx >= n) { idx = 0; }
        }
      for (idx = c;;)
        {
        poly1.push_back(innerPoly[idx]);
        if (idx == d) { break; }
        if (++idx >= m) { idx = 0; }
        }

      // Generate poly2
      vtkClipOutlinePoly poly2;
      for (idx = b;;)
        {
        poly2.push_back(outerPoly[idx]);
        if (idx == a) { break; }
        if (++idx >= n) { idx = 0; }
        }
      for (idx = d;;)
        {
        poly2.push_back(innerPoly[idx]);
        if (idx == c) { break; }
        if (++idx >= m) { idx = 0; }
        }

      // Replace outerPoly and innerPoly with these new polys
      polys[outerPolyId] = poly1;
      polys[innerPolyId] = poly2;

      // Move innerPolyId into its own group
      polyGroup.erase(polyGroup.begin()+1);
      polyGroups[innerPolyId].push_back(innerPolyId);

      // If there are other interior polys in the group, find out whether
      // they are in poly1 or poly2
      if (polyGroup.size() > 1)
        {
        double *pp = new double[3*poly1.size()];
        double bounds[6];
        double tol2;
        vtkClipOutlinePrepareForPolyInPoly(poly1, points, pp, bounds, tol2);

        size_t ii = 1;
        while (ii < polyGroup.size())
          {
          if (vtkClipOutlinePolyInPoly(poly1, polys[polyGroup[ii]],
                                       points, normal, pp, bounds, tol2))
            {
            // Keep this poly in polyGroup
            ii++;
            }
          else
            {
            // Move this poly to poly2 group
            polyGroups[innerPolyId].push_back(polyGroup[ii]);
            polyGroup.erase(polyGroup.begin()+ii);

            // Reduce the groupId to ensure that this new group
            // will get cut
            if (innerPolyId < groupId)
              {
              groupId = innerPolyId;
              }
            }
          }
        delete [] pp;

        // Continue without incrementing groupId
        continue;
        }
      }

    // Increment to the next group
    groupId++;
    }
}
