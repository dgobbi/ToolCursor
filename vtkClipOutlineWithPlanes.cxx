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

vtkCxxRevisionMacro(vtkClipOutlineWithPlanes, "$Revision: 1.3 $");
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

  // Skip through the old lines
  vtkIdType npts, *pts;
  lines->InitTraversal();
  for (vtkIdType cellId = 0; cellId < firstLine; cellId++)
    {
    lines->GetNextCell(npts, pts);
    }

  // Save the location of the first new line
  vtkIdType saveLoc = lines->GetTraversalLocation();

  // Temporary storage for all the new polys
  vtkstd::vector<vtkstd::vector<vtkIdType> > newPolys;

  // ---------------------------------------------------
  // Here is the code for creating polygons from line segments

  // Bitfield for marking lines as used
  vtkClipOutlineBitArray usedLines;

  size_t numNewPolys = 0;
  vtkIdType remainingLines = numNewLines;
  while (remainingLines > 0)
    {
    // Create a new poly
    newPolys.resize(++numNewPolys);
    vtkstd::vector<vtkIdType> &poly = newPolys[numNewPolys-1];

    int completePoly = 0;
    int noLinesMatch = 0;
    while (!completePoly && !noLinesMatch && remainingLines > 0)
      {
      noLinesMatch = 1;
      lines->SetTraversalLocation(saveLoc);
      for (vtkIdType lineId = 0; lineId < numNewLines; lineId++)
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

  // After the above loop, all the lines have been collected
  // into connected groups.  If we are lucky, those
  // connected groups are simple, convex polygons.
  // However, there are other possibilities:
  
  // Some groups might be self-intersecting.  This is easy
  // to detect:  if the same point appears multiple times,
  // the group must be split at that point.

  // Some of these groups may be "holes" in a larger,
  // enclosing polygon.  For each polygon, these holes
  // can be found by searching for other polygons within
  // the enclosing polygon.  Then the interior polygon can be
  // joined with the enclosing polygon by cutting a keyhole.
  // The keyhole should go between the closest points on the
  // interior and exterior polygon, but must not intersect
  // any line segments.  The sense of the internal polygon
  // will have to be reversed.

  // Some of these groups might be concave polygons.
  // All polygons must be converted into quads and triangles
  // to guarantee that all polys are convex.

  // ---------------------------------------------------
  // Check for self-intersection. Split the figure-eights.

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

  // ---------------------------------------------------
  // Correct the sense of the polygons

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

  // ---------------------------------------------------
  // Check for polygons within polygons.  Group the polygons
  // if they are within each other.  Reverse the sense of 
  // the interior "hole" polygons.  A hole within a hole
  // will be reversed twice.

  // Make an array of bitsets.  Start with one poly per bitset,
  // but consolidate when polys are found inside other polys.
  vtkstd::vector<vtkClipOutlineBitArray> polyGroups;
  polyGroups.resize(numNewPolys);
  for (size_t i = 0; i < numNewPolys; i++)
    {
    polyGroups[i].set(i, 1);
    }

/*
  for (size_t i = 0; i < numNewPolys; i++)
    {
    size_t n = newPolys[i].size();

    if (n < 3) { continue; }

    double *pp = new double[3*n];
    double bounds[6];
    bounds[0] = bounds[1] = bounds[2] = 0;
    bounds[3] = bounds[4] = bounds[5] = 0;

    for (size_t k = 0; k < n; k++)
      {
      double *p = &pp[3*k];
      points->GetPoint(newPolys[i][k], p);
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
    double tol2 = (ps[0]*ps[0] + ps[1]*ps[1] + ps[2]*ps[2])*(1e-5 * 1e-5);

    for (size_t j = 0; j < numNewPolys; j++)
      {
      size_t m = newPolys[j].size();
      if (j == i || m < 3) { continue; }

      // Find a vertex of poly "j" that isn't on the edge of poly "i".
      // This is necessary or the PointInPolygon might return "true"
      // based only on roundoff error.

      double p[3];
      int allPointsOnEdges = 1;

      for (size_t jj = 0; jj < m; jj++)
        {          
        points->GetPoint(newPolys[j][jj], p);

        int pointOnEdge = 0;
        double q1[3], q2[3];
        points->GetPoint(newPolys[i][n-1], q1);
        for (size_t ii = 0; ii < n; ii++)
          {
          points->GetPoint(newPolys[i][ii], q2);
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

      // There should also be a check to see if all the verts match.
      // If they do, both polys should be removed.

      if (allPointsOnEdges || 
          vtkPolygon::PointInPolygon(p, n, pp, bounds, normal) == 1)
        {
        // Reverse the interior polygon, and mark it as being grouped
        // with the outer polygon.  Together, grouped polygons can
        // be used as constraints for a Delaunay triangulation.

        size_t m2 = m/2;
        for (size_t kk = 0; kk < m2; kk++)
          {
          vtkIdType tmpId = newPolys[j][kk];
          newPolys[j][kk] = newPolys[j][m - kk - 1];
          newPolys[j][m - kk - 1] = tmpId;
          }

        // Search through existing polygon groups (stored as bitsets)
        int foundGroup = 0;
        size_t firstGroupId = 0;

        for (size_t groupId = 0; groupId < polyGroups.size(); groupId++)
          {
          if (polyGroups[groupId].get(i) || polyGroups[groupId].get(j))
            {
            if (!foundGroup)
              {
              foundGroup = 1;
              firstGroupId = groupId;
              polyGroups[firstGroupId].set(i, 1);
              polyGroups[firstGroupId].set(j, 1);
              }
            else
              {
              polyGroups[firstGroupId].merge(polyGroups[groupId]);
              polyGroups[groupId] = polyGroups[polyGroups.size()-1];
              polyGroups.resize(polyGroups.size()-1);
              groupId--;
              }
            } 
          } 
        }
      }
      delete [] pp;
    }
*/
  // Now that the inside-outside check has been done, almost
  // everything is ready for a Delaunay triangulation.

  // The Delaunay requires creating a new set of points that
  // are transformed to the 2D plane, and that is limited only
  // to the points that form the polygons to be triangulated.
  // The clipping plane normal and position will be used to
  // compute the transform.  A vector that maps the new pointIds
  // to the original pointIds will be needed in order to
  // merge the triangulation with our output cell array.

  // Delaunay furthermore requires a data set with a group's
  // worth of polygons, with all the Ids converted to the new Ids.
  // This polydata data set will be used for both of the inputs
  // to Delaunay2D, i.e. the points and the contraints.

  // The output of Delaunay will be the new polys.  The pointIds
  // must be converted back when they are stored on our polys list.

  // ---------------------------------------------------
  // Add the polygons to the vtkCellArray.  Triangulate where necessary.

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
    // Make a map for recovering pointIds after triangulation has finished
    vtkstd::vector<vtkIdType> pointMap;
    vtkIdType newPointId = 0;

    sourcePolys->Initialize();
    sourcePoints->Initialize();
    sourceData->Initialize();
    sourceData->SetPoints(sourcePoints);
    sourceData->SetPolys(sourcePolys);

    // Go through all polys, check whether they are in the group.
    for (size_t i = 0; i < numNewPolys; i++)
      {
      size_t n = newPolys[i].size();

      // If the poly is not in the group or is a line, then skip it
      if (polyGroups[groupId].get(i) == 0 || n < 3)
        {
        continue;
        }

      // Create the source polys for a constrained Delaunay
      sourcePolys->InsertNextCell(n);
      for (size_t j = 0; j < n; j++)
        {
        // Use homogeneous point with Matrix4x4
        double point[4];
        point[3] = 1.0;

        vtkIdType pointId = newPolys[i][j];
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


