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

vtkCxxRevisionMacro(vtkClipOutlineWithPlanes, "$Revision: 1.1 $");
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
}

//----------------------------------------------------------------------------
vtkClipOutlineWithPlanes::~vtkClipOutlineWithPlanes()
{
  if (this->ClippingPlanes)
    {
    this->ClippingPlanes->Delete();
    this->ClippingPlanes = 0;
    }
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

  // The points.
  vtkPoints *points = vtkPoints::New();
  points->DeepCopy(input->GetPoints());

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
    if (lineScalars)
      {
      lineScalars->Delete();
      }
    if (polys)
      {
      polys->Delete();
      }
    if (polyScalars)
      {
      polyScalars->Delete();
      }

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
  vtkIncrementalPointLocator *locator = vtkMergePoints::New();
  vtkPoints *newPoints = vtkPoints::New();

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
    this->ClipCells(points, pointScalars, locator, 1, lines, newLines,
                    inPointData, outPointData, inLineData, outLineData);

    if (polys)
      {
      // Get the number of lines remaining after the clipping
      vtkIdType numClipLines = newLines->GetNumberOfCells();

      // Cut the polys to generate more lines
      this->ContourCells(points, pointScalars, locator, 2, polys, newLines,
                         inPointData, outPointData, inPolyData, outLineData);
      
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

      // Clip the polys
      this->ClipCells(points, pointScalars, locator, 2, polys, newPolys,
                      inPointData, outPointData, inPolyData, outPolyData);

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

  newPoints->Delete();
  locator->Delete();
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
void vtkClipOutlineWithPlanes::ClipCells(
  vtkPoints *points, vtkDoubleArray *pointScalars,
  vtkIncrementalPointLocator *locator, int dimensionality,
  vtkCellArray *inputCells, vtkCellArray *outputCells,
  vtkPointData *inPD, vtkPointData *outPD,
  vtkCellData *inCD, vtkCellData *outCD)
{
  vtkDoubleArray *cellPointScalars = vtkDoubleArray::New();
  vtkGenericCell *cell = vtkGenericCell::New();

  vtkIdType numCells = inputCells->GetNumberOfCells();
  inputCells->InitTraversal();
  for (vtkIdType cellId = 0; cellId < numCells; cellId++)
    {
    vtkIdType numPts, *pts;
    inputCells->GetNextCell(numPts, pts);

    // Set the cell type from the dimensionality
    if (dimensionality == 1)
      {
      if (numPts == 2) { cell->SetCellTypeToLine(); }
      else { cell->SetCellTypeToPolyLine(); }
      }
    else // (dimensionality == 2)
      {
      if (numPts == 3) { cell->SetCellTypeToTriangle(); }
      else if (numPts == 4) { cell->SetCellTypeToQuad(); }
      else { cell->SetCellTypeToPolygon(); }
      }

    vtkPoints *cellPts = cell->GetPoints();
    vtkIdList *cellIds = cell->GetPointIds();

    cellPts->SetNumberOfPoints(numPts);
    cellIds->SetNumberOfIds(numPts);
    cellPointScalars->SetNumberOfValues(numPts);

    // Copy everything over to the temporary cell
    for (vtkIdType i = 0; i < numPts; i++)
      {
      double point[3];
      points->GetPoint(pts[i], point);
      cellPts->SetPoint(i, point);
      cellIds->SetId(i, pts[i]);
      double s = pointScalars->GetValue(cellIds->GetId(i));
      cellPointScalars->SetValue(i, s);
      }

    cell->Clip(0, cellPointScalars, locator, outputCells, inPD, outPD,
               inCD, cellId, outCD, 0);
    }

  cellPointScalars->Delete();
  cell->Delete();
}             

//----------------------------------------------------------------------------
void vtkClipOutlineWithPlanes::ContourCells(
  vtkPoints *points, vtkDoubleArray *pointScalars,
  vtkIncrementalPointLocator *locator, int dimensionality,
  vtkCellArray *inputCells, vtkCellArray *outputCells,
  vtkPointData *inPD, vtkPointData *outPD,
  vtkCellData *inCD, vtkCellData *outCD)
{
  vtkDoubleArray *cellPointScalars = vtkDoubleArray::New();
  vtkGenericCell *cell = vtkGenericCell::New();
  vtkCellArray *verts = vtkCellArray::New();

  vtkIdType numCells = inputCells->GetNumberOfCells();
  inputCells->InitTraversal();
  for (vtkIdType cellId = 0; cellId < numCells; cellId++)
    {
    vtkIdType numPts, *pts;
    inputCells->GetNextCell(numPts, pts);

    // Set the cell type from the dimensionality
    if (dimensionality == 1)
      {
      if (numPts == 2) { cell->SetCellTypeToLine(); }
      else { cell->SetCellTypeToPolyLine(); }
      }
    else // (dimensionality == 2)
      {
      if (numPts == 3) { cell->SetCellTypeToTriangle(); }
      else if (numPts == 4) { cell->SetCellTypeToQuad(); }
      else { cell->SetCellTypeToPolygon(); }
      }

    vtkPoints *cellPts = cell->GetPoints();
    vtkIdList *cellIds = cell->GetPointIds();

    cellPts->SetNumberOfPoints(numPts);
    cellIds->SetNumberOfIds(numPts);
    cellPointScalars->SetNumberOfValues(numPts);

    // Copy everything over to the temporary cell
    for (vtkIdType i = 0; i < numPts; i++)
      {
      double point[3];
      points->GetPoint(pts[i], point);
      cellPts->SetPoint(i, point);
      cellIds->SetId(i, pts[i]);
      double s = pointScalars->GetValue(cellIds->GetId(i));
      cellPointScalars->SetValue(i, s);
      }

    cell->Contour(0, cellPointScalars, locator, verts, outputCells, 0,
                  inPD, outPD, inCD, cellId, outCD);
    }

  verts->Delete();
  cellPointScalars->Delete();
  cell->Delete();
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
// A helper class: a bitfield that is always as large as needed

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

  for (size_t i = 0; i < numNewPolys; i++)
    {
    size_t n = newPolys[i].size();

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

    // Triangularize any polys with more than 4 points
    if (n > 4)
      {
      vtkPolygon *polygon = vtkPolygon::New();
      vtkIdList *triangles = vtkIdList::New();

      polygon->Points->SetNumberOfPoints(n);
      polygon->PointIds->SetNumberOfIds(n);

      for (size_t j = 0; j < n; j++)
        {
        vtkIdType pointId = newPolys[i][j];
        double point[3];
        points->GetPoint(pointId, point);
        polygon->Points->SetPoint(static_cast<vtkIdType>(j), point);
        polygon->PointIds->SetId(static_cast<vtkIdType>(j), pointId);
        }

      polygon->Triangulate(triangles);

      vtkUnsignedCharArray *scalars =
        vtkUnsignedCharArray::SafeDownCast(outCD->GetScalars());

      vtkIdType m = triangles->GetNumberOfIds();
      for (vtkIdType k = 0; k < m; k += 3)
        {
        vtkIdType cellId = polys->InsertNextCell(3);
        polys->InsertCellPoint(newPolys[i][triangles->GetId(k + 0)]);
        polys->InsertCellPoint(newPolys[i][triangles->GetId(k + 1)]);
        polys->InsertCellPoint(newPolys[i][triangles->GetId(k + 2)]);

        if (scalars)
          {
          scalars->InsertTupleValue(cellId, color);
          }
        }

      polygon->Delete();
      triangles->Delete();
      }
    else if (n > 2)
      {
      // Add any triangles or quads directly

      vtkIdType cellId = polys->InsertNextCell(static_cast<vtkIdType>(n));

      for (size_t j = 0; j < n; j++)
        {
        polys->InsertCellPoint(newPolys[i][j]);
        }

      vtkUnsignedCharArray *scalars =
        vtkUnsignedCharArray::SafeDownCast(outCD->GetScalars());

      if (scalars)
        {
        scalars->InsertTupleValue(cellId, color);
        }
      }
    }
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


