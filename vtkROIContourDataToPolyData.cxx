/*=========================================================================

  Program:   ToolCursor
  Module:    vtkROIContourDataToPolyData.cxx

  Copyright (c) 2010 David Gobbi
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkROIContourDataToPolyData.h"

#include "vtkROIContourData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkImageData.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkIntArray.h"
#include "vtkDoubleArray.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkPlane.h"
#include "vtkKochanekSpline.h"
#include "vtkMath.h"

vtkStandardNewMacro(vtkROIContourDataToPolyData);
vtkCxxSetObjectMacro(vtkROIContourDataToPolyData,SelectionPlane,vtkPlane);
vtkCxxSetObjectMacro(vtkROIContourDataToPolyData,Spline,vtkSpline);

//----------------------------------------------------------------------------
vtkROIContourDataToPolyData::vtkROIContourDataToPolyData()
{
  this->SelectionPlane = NULL;
  this->SelectionPlaneTolerance = 0.5;
  this->Subdivision = 0;
  this->SubdivisionTarget = 1.0;
  this->Spline = 0;

  this->SplineX = 0;
  this->SplineY = 0;
  this->SplineZ = 0;
  this->KnotPositions = 0;
}

//----------------------------------------------------------------------------
vtkROIContourDataToPolyData::~vtkROIContourDataToPolyData()
{
  if (this->SelectionPlane)
    {
    this->SelectionPlane->Delete();
    }
  if (this->Spline)
    {
    this->Spline->Delete();
    }
  if (this->SplineX)
    {
    this->SplineX->Delete();
    this->SplineY->Delete();
    this->SplineZ->Delete();
    }
  if (this->KnotPositions)
    {
    this->KnotPositions->Delete();
    }
}

//----------------------------------------------------------------------------
void vtkROIContourDataToPolyData::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "SelectionPlane: " << this->SelectionPlane << "\n";
  os << indent << "SelectionPlaneTolerance: "
     << this->SelectionPlaneTolerance << "\n";
  os << indent << "Subdivision: " << (this->Subdivision ? "On\n" : "Off\n");
  os << indent << "SubdivisionTarget: " << this->SubdivisionTarget << "\n";
  os << indent << "Spline: " << this->Spline << "\n";
  
}

//----------------------------------------------------------------------------
int vtkROIContourDataToPolyData::FillInputPortInformation(
  int, vtkInformation *info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkROIContourData");
  return 1;
}

//----------------------------------------------------------------------------
int vtkROIContourDataToPolyData::ComputePipelineMTime(
  vtkInformation* vtkNotUsed(request),
  vtkInformationVector** vtkNotUsed(inputVector),
  vtkInformationVector* vtkNotUsed(outputVector),
  int vtkNotUsed(requestFromOutputPort),
  unsigned long* mtime)
{
  unsigned long mTime = this->GetMTime();

  vtkPlane *plane = this->SelectionPlane;

  if (plane)
    {
    unsigned long planeMTime = plane->GetMTime();
    if (planeMTime > mTime)
      {
      mTime = planeMTime;
      }
    }

  *mtime = mTime;

  return 1;
}

//----------------------------------------------------------------------------
// Compute SplineX, SplineY, and SplineZ for the given points, using the
// accumulated distance between points as the parameter for the spline.
// Set "closed" to one if the splines should form a closed loop.
// The reference parameter tmax will be set to the maximum parameter
// value for the spline, and dmax will be set to the full length of
// the contour, computed by summing the lengths of the line segments.
void vtkROIContourDataToPolyData::ComputeSpline(
  vtkPoints *points, int closed, double &tmax, double &dmax)
{
  if (this->Spline && this->SplineX &&
      !this->SplineX->IsA(this->Spline->GetClassName()))
    {
    this->SplineX->Delete();
    this->SplineY->Delete();
    this->SplineZ->Delete();

    this->SplineX = 0;
    this->SplineY = 0;
    this->SplineZ = 0;
    }
  if (this->SplineX == 0)
    {
    if (this->Spline)
      {
      this->SplineX = this->Spline->NewInstance();
      this->SplineX->DeepCopy(this->Spline);
      this->SplineY = this->Spline->NewInstance();
      this->SplineY->DeepCopy(this->Spline);
      this->SplineZ = this->Spline->NewInstance();
      this->SplineZ->DeepCopy(this->Spline);
      }
    else
      {
      this->SplineX = vtkKochanekSpline::New();
      this->SplineY = vtkKochanekSpline::New();
      this->SplineZ = vtkKochanekSpline::New();
      }
    }
  if (this->KnotPositions == 0)
    {
    this->KnotPositions = vtkDoubleArray::New();
    }

  vtkSpline *xspline = this->SplineX;
  vtkSpline *yspline = this->SplineY;
  vtkSpline *zspline = this->SplineZ;
  vtkDoubleArray *knots = this->KnotPositions;

  // initialize the spline
  xspline->RemoveAllPoints();
  yspline->RemoveAllPoints();
  zspline->RemoveAllPoints();
  knots->Initialize();

  // set whether splines are closed
  xspline->SetClosed(closed);
  yspline->SetClosed(closed);
  zspline->SetClosed(closed);

  // get the number of points
  vtkIdType n = points->GetNumberOfPoints();
  double p0[3], p[3];

  // factor between real distance and parametric distance
  double f = 1.0;
  // the length of the implicit segment for closed loops
  double lastd = 0;

  // verify that there are enough knots for the spline
  if (n < 2)
    {
    tmax = 0;
    dmax = 0;
    return;
    }

  // get the first and last point
  points->GetPoint(0, p0);

  if (closed)
    {
    // require a tolerance, base it off the desired subdivision
    double tol = this->SubdivisionTarget*1e-3;
    tol *= tol;

    // ignore last point if same as first point
    vtkIdType m = n;
    do
      {
      points->GetPoint(--m, p);
      lastd = vtkMath::Distance2BetweenPoints(p0, p);
      }
    while (m > 0 && lastd < tol);
    n = m + 1;

    // set factor to scale the implicit segment to unity
    if (lastd > 0)
      {
      lastd = sqrt(lastd);
      f = 1.0/lastd;
      }
    }

  // verify that there are still enough knots for the spline
  if (n < 2)
    {
    tmax = 0;
    dmax = 0;
    return;
    }

  // add all the points to the spline
  double d = 0.0;
  for (vtkIdType i = 0; i < n; i++)
    {
    points->GetPoint(i, p);
    d += sqrt(vtkMath::Distance2BetweenPoints(p0, p));

    double t = f*d;

    xspline->AddPoint(t, p[0]);
    yspline->AddPoint(t, p[1]);
    zspline->AddPoint(t, p[2]);
    knots->InsertNextValue(t);

    p0[0] = p[0];
    p0[1] = p[1];
    p0[2] = p[2];
    }

  // do the spline precomputations
  xspline->Compute();
  yspline->Compute();
  zspline->Compute();

  // the spline is valid over t = [0, tmax]
  d += lastd;
  tmax = f*d;
  dmax = d;

  // add another knot point for closed splines
  if (closed)
    {
    knots->InsertNextValue(tmax);
    }
}


//----------------------------------------------------------------------------
void vtkROIContourDataToPolyData::GenerateSpline(
  vtkPoints *contourPoints, int closed,
  vtkPoints *points, vtkCellArray *lines, vtkIntArray *subIds)
{
  double tmax, dmax;
  this->ComputeSpline(contourPoints, closed, tmax, dmax);

  vtkDoubleArray *knots = this->KnotPositions;
  vtkIdType m = knots->GetNumberOfTuples();

  if (m > 1)
    {
    vtkSpline *xspline = this->SplineX;
    vtkSpline *yspline = this->SplineY;
    vtkSpline *zspline = this->SplineZ;
    vtkIdType id0 = points->GetNumberOfPoints();
    double t0 = 0;
    double f = dmax/(tmax*this->SubdivisionTarget);
    for (vtkIdType j = 1; j < m; j++)
      {
      double t1 = knots->GetValue(j);
      int n = vtkMath::Floor((t1 - t0)*f) + 1;
      for (int i = 0; i < n; i++)
        {
        double t = (t0*(n-i) + t1*i)/n;
        double p[3];
        p[0] = xspline->Evaluate(t);
        p[1] = yspline->Evaluate(t);
        p[2] = zspline->Evaluate(t);
        points->InsertNextPoint(p);
        if (subIds)
          {
          subIds->InsertNextValue(static_cast<int>(j-1));
          }
        }
      t0 = t1;
      }

    if (!closed)
      {
      double p[3];
      p[0] = xspline->Evaluate(tmax);
      p[1] = yspline->Evaluate(tmax);
      p[2] = zspline->Evaluate(tmax);
      points->InsertNextPoint(p);
      if (subIds)
        {
        subIds->InsertNextValue(m-1);
        }
      }

    vtkIdType id1 = points->GetNumberOfPoints();
    lines->InsertNextCell(id1 - id0 + (closed != 0));
    for (vtkIdType id = id0; id < id1; id++)
      {
      lines->InsertCellPoint(id);
      }
    if (closed)
      {
      lines->InsertCellPoint(id0);
      }
    }
  knots->Initialize();
}

//----------------------------------------------------------------------------
int vtkROIContourDataToPolyData::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // Get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // Get the input and output
  vtkROIContourData *input = vtkROIContourData::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData *output = vtkPolyData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  // For selecting which contours to include in the output
  vtkPlane *plane = this->SelectionPlane;
  double tol = this->SelectionPlaneTolerance;

  // The output points and cells
  vtkPoints *outPoints = vtkPoints::New(VTK_DOUBLE);
  vtkCellArray *lines = 0;
  vtkCellArray *verts = 0;

  // The scalars
  vtkIntArray *contourIds = vtkIntArray::New();
  contourIds->SetName("Labels");
  vtkIntArray *contourSubIds = vtkIntArray::New();
  contourSubIds->SetName("SubIds");

  // Go through all the contours
  int n = input->GetNumberOfContours();
  for (int i = 0; i < n; i++)
    {
    vtkPoints *points = input->GetContourPoints(i);
    int t = input->GetContourType(i);

    contourIds->InsertNextValue(i);

    if (points)
      {
      vtkIdType m = points->GetNumberOfPoints();
      bool includeContour = true;

      if (plane)
        {
        for (int j = 0; j < m; j++)
          {
          double p[3];
          points->GetPoint(j, p);
          double d = plane->DistanceToPlane(p);
          if (d < -tol || d > tol)
            {
            includeContour = false;
            break;
            }
          }
        }

      // Include this contour in the output
      if (includeContour)
        {
        int closed = 0;
        vtkCellArray *cells = 0;
        vtkIdType cellSize = m;
        if (t == vtkROIContourData::POINT)
          {
          if (!verts)
            {
            verts = vtkCellArray::New();
            }
          cells = verts;
          }
        else
          {
          if (!lines)
            {
            lines = vtkCellArray::New();
            }
          cells = lines;

          if (t == vtkROIContourData::CLOSED_PLANAR)
            {
            cellSize = m + 1;
            closed = 1;
            }
          }

        if (this->Subdivision && m > 2 && t != vtkROIContourData::POINT)
          {
          // Use a spline to subdivide and smooth contour
          this->GenerateSpline(
            points, closed, outPoints, lines, contourSubIds);
          }
        else
          {
          // Add the contour without subdivision
          cells->InsertNextCell(cellSize);
          vtkIdType firstPointId = outPoints->GetNumberOfPoints();

          for (int j = 0; j < m; j++)
            {
            double p[3];
            points->GetPoint(j, p);
            vtkIdType pointId = outPoints->InsertNextPoint(p);
            if (contourSubIds)
              {
              contourSubIds->InsertNextValue(j);
              }
            cells->InsertCellPoint(pointId);
            }

          // Close the contour, if necessary
          if (cellSize > m)
            {
            cells->InsertCellPoint(firstPointId);
            }
          }
        }
      }
    }


  output->SetPoints(outPoints);
  output->SetLines(lines);
  output->SetVerts(verts);
  output->GetCellData()->SetScalars(contourIds);
  output->GetPointData()->AddArray(contourSubIds);

  if (outPoints)
    {
    outPoints->Delete();
    }
  if (lines)
    {
    lines->Delete();
    }
  if (verts)
    {
    verts->Delete();
    }
  if (contourIds)
    {
    contourIds->Delete();
    }
  if (contourSubIds)
    {
    contourSubIds->Delete();
    }

  // assign colors to the output points
  unsigned char color[3] = { 255, 0, 0 };
  vtkUnsignedCharArray *colors = vtkUnsignedCharArray::New();
  colors->SetNumberOfComponents(3);
  colors->SetName("Colors");
  vtkIdType m = outPoints->GetNumberOfPoints();
  colors->SetNumberOfTuples(m);
  for (vtkIdType j = 0; j < m; j++)
    {
    colors->SetTupleValue(j, color);
    }

  output->GetPointData()->SetScalars(colors);
  colors->Delete();

  return 1;
}
