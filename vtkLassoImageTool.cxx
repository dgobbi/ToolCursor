/*=========================================================================

  Program:   ToolCursor
  Module:    vtkLassoImageTool.cxx

  Copyright (c) 2010 David Gobbi
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkLassoImageTool.h"
#include "vtkObjectFactory.h"

#include "vtkROIContourData.h"
#include "vtkROIContourDataToPolyDataFilter.h"
#include "vtkGeometricCursorShapes.h"
#include "vtkToolCursor.h"
#include "vtkCamera.h"
#include "vtkRenderer.h"
#include "vtkMatrix4x4.h"
#include "vtkMath.h"
#include "vtkCardinalSpline.h"
#include "vtkGlyph3D.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkActor.h"
#include "vtkDataSetMapper.h"
#include "vtkProperty.h"

#include "vtkVolumePicker.h"

vtkStandardNewMacro(vtkLassoImageTool);

//----------------------------------------------------------------------------
vtkLassoImageTool::vtkLassoImageTool()
{
  this->Matrix = vtkMatrix4x4::New();

  vtkGeometricCursorShapes *shapes = vtkGeometricCursorShapes::New();

  this->SplineX = vtkCardinalSpline::New();
  this->SplineY = vtkCardinalSpline::New();

  this->ROIData = vtkROIContourData::New();

  this->ROIDataToPolyData = vtkROIContourDataToPolyDataFilter::New();
  this->ROIDataToPolyData->SetInput(this->ROIData);

  this->Glyph3D = vtkGlyph3D::New();
  this->Glyph3D->SetColorModeToColorByScalar();
  this->Glyph3D->SetScaleFactor(0.7);
  this->Glyph3D->SetInputConnection(this->ROIDataToPolyData->GetOutputPort());
  this->Glyph3D->SetSource(vtkPolyData::SafeDownCast(
    shapes->GetShapeData("Sphere")));

  this->GlyphMapper = vtkDataSetMapper::New();
  this->GlyphMapper->SetInputConnection(this->Glyph3D->GetOutputPort());

  this->GlyphActor = vtkActor::New();
  this->GlyphActor->PickableOff();
  this->GlyphActor->SetMapper(this->GlyphMapper);
  this->GlyphActor->GetProperty()->SetColor(1,0,0);
  this->GlyphActor->SetUserMatrix(this->Matrix);

  this->ContourData = vtkPolyData::New();
  vtkPoints *linePoints = vtkPoints::New();
  vtkCellArray *lines = vtkCellArray::New();
  this->ContourData->SetPoints(linePoints);
  this->ContourData->SetLines(lines);
  linePoints->Delete();
  lines->Delete();

  this->ContourMapper = vtkDataSetMapper::New();
  this->ContourMapper->SetInput(this->ContourData);

  this->ContourActor = vtkActor::New();
  this->ContourActor->PickableOff();
  this->ContourActor->SetMapper(this->ContourMapper);
  this->ContourActor->GetProperty()->SetColor(1,0,0);
  this->ContourActor->SetUserMatrix(this->Matrix);

  shapes->Delete();

  this->CurrentPointId = -1;
  this->InitialPointPosition[0] = 0;
  this->InitialPointPosition[1] = 1;
  this->InitialPointPosition[2] = 2;
}

//----------------------------------------------------------------------------
vtkLassoImageTool::~vtkLassoImageTool()
{
  if (this->ROIData)
    {
    this->ROIData->Delete();
    }

  this->SplineX->Delete();
  this->SplineY->Delete();

  this->Matrix->Delete();
  this->ROIDataToPolyData->Delete();
  this->ContourData->Delete();
  this->Glyph3D->Delete();
  this->GlyphActor->Delete();
  this->GlyphMapper->Delete();
}

//----------------------------------------------------------------------------
void vtkLassoImageTool::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

//----------------------------------------------------------------------------
void vtkLassoImageTool::SetMarker(vtkPolyData *data)
{
  if (data != this->Glyph3D->GetSource())
    {
    this->Glyph3D->SetSource(data);
    this->Modified();
    }
}

//----------------------------------------------------------------------------
vtkPolyData *vtkLassoImageTool::GetMarker()
{
  return this->Glyph3D->GetSource();
}

//----------------------------------------------------------------------------
void vtkLassoImageTool::StartAction()
{
  this->Superclass::StartAction();

  vtkToolCursor *cursor = this->GetToolCursor();

  double position[3];
  cursor->GetPosition(position);

  vtkROIContourData *data = this->ROIData;
  vtkPoints *points = 0;
  if (data->GetNumberOfContours() == 0)
    {
    points = vtkPoints::New();
    this->ROIData->SetNumberOfContours(1);
    this->ROIData->SetContourPoints(0, points);
    points->Delete();
    }
  else
    {
    points= this->ROIData->GetContourPoints(0);
    }

  double tol = 3;
  // Is the mouse over a previous point?
  this->CurrentPointId = -1;
  vtkIdType n = points->GetNumberOfPoints();
  for (vtkIdType i = 0; i < n; i++)
    {
    double p[3];
    points->GetPoint(i, p);
    if (vtkMath::Distance2BetweenPoints(position, p) < tol*tol)
      {
      this->CurrentPointId = i;
      this->InitialPointPosition[0] = p[0];
      this->InitialPointPosition[1] = p[1];
      this->InitialPointPosition[2] = p[2];
      break;
      }
    }

  if (this->CurrentPointId < 0)
    {
    points->InsertNextPoint(position);
    this->ROIData->Modified();
    this->UpdateContourData();
    cursor->GetRenderer()->AddViewProp(this->GlyphActor);
    cursor->GetRenderer()->AddViewProp(this->ContourActor);

    cerr << "inserting point " << points->GetNumberOfPoints() << " " << position[0] << ", " << position[1] << ", " << position[2] << "\n";
    }
}

//----------------------------------------------------------------------------
void vtkLassoImageTool::StopAction()
{
  this->Superclass::StopAction();
}

//----------------------------------------------------------------------------
void vtkLassoImageTool::DoAction()
{
  this->Superclass::DoAction();

  vtkToolCursor *cursor = this->GetToolCursor();

  // Get the initial point.
  double p0[3];
  this->GetStartPosition(p0);

  // Get the current position
  double position[3];
  cursor->GetPosition(position);

  double dx = position[0] - p0[0];
  double dy = position[1] - p0[1];

  vtkPoints *points = 0;
  if (this->ROIData->GetNumberOfContours() > 0)
    {
    points = this->ROIData->GetContourPoints(0);
    }

  if (points && this->CurrentPointId >= 0)
    {
    double p[3];
    p[0] = this->InitialPointPosition[0] + dx;
    p[1] = this->InitialPointPosition[1] + dy;
    p[2] = this->InitialPointPosition[2];
    points->SetPoint(this->CurrentPointId, p);

    this->ROIData->Modified();
    this->UpdateContourData();
    }
}

//----------------------------------------------------------------------------
void vtkLassoImageTool::UpdateSpline(double &tmax, double &dmax)
{

  vtkSpline *xspline = this->SplineX;
  vtkSpline *yspline = this->SplineY;

  // initialize the spline
  xspline->RemoveAllPoints();
  yspline->RemoveAllPoints();
  xspline->ClosedOff();
  yspline->ClosedOff();

  // get the number of points and the first/last point
  vtkPoints *points = 0;
  if (this->ROIData->GetNumberOfContours() > 0)
    {
    points = this->ROIData->GetContourPoints(0);
    }

  vtkIdType n = points->GetNumberOfPoints();
  double p[3];
  double p0[2], p1[2];

  int xj = 0;
  int yj = 1;

  points->GetPoint(n-1, p);
  p0[0] = p[xj];
  p0[1] = p[yj];

  points->GetPoint(0, p);
  p1[0] = p[xj];
  p1[1] = p[yj];

  // factor between real distance and parametric distance
  double f = 1.0;
  // the length of the implicit segment for closed loops
  double lastd = 0;

  // if first and last point are same, spline is closed
  double dx = p1[0] - p0[0];
  double dy = p1[1] - p0[1];
  double d2 = dx*dx + dy*dy;
  cerr << "d2 " << d2 << "\n";
  while (d2 < 100 && n > 1)
    {
    cerr << "closing\n";
    n -= 1;
    points->GetPoint(n-1, p);
    p0[0] = p[xj];
    p0[1] = p[yj];

    xspline->ClosedOn();
    yspline->ClosedOn();

    // vtkSpline considers the parametric length of the implicit
    // segment of closed loops to be unity, so set "f" so that
    // multiplying the real length of that segment by "f" gives unity.
    dx = p1[0] - p0[0];
    dy = p1[1] - p0[1];
    d2 = dx*dx + dy*dy;
    lastd = sqrt(d2);
    if (lastd > 0)
      {
      f = 1.0/lastd;
      }
    }

  // Add all the points to the spline.
  double d = 0.0;
  for (vtkIdType i = 0; i < n; i++)
    {
    p0[0] = p1[0];
    p0[1] = p1[1];

    points->GetPoint(i, p);
    p1[0] = p[xj];
    p1[1] = p[yj];

    dx = p1[0] - p0[0];
    dy = p1[1] - p0[1];

    d += sqrt(dx*dx + dy*dy);

    double t = f*d;

    xspline->AddPoint(t, p1[0]);
    yspline->AddPoint(t, p1[1]);
    }

  // Do the spline precomputations
  xspline->Compute();
  yspline->Compute();

  // The spline is valid over t = [0, tmax]
  d += lastd;
  tmax = f*d;
  dmax = d;
}

//----------------------------------------------------------------------------
void vtkLassoImageTool::UpdateContourData()
{
  vtkPoints *points = this->ContourData->GetPoints();
  vtkCellArray *lines = this->ContourData->GetLines();

  vtkPoints *knots = 0;
  if (this->ROIData->GetNumberOfContours() > 0)
    {
    knots = this->ROIData->GetContourPoints(0);
    }

  if (!knots)
    {
    return;
    }

  double z = 0;
  if (knots->GetNumberOfPoints() > 0)
    {
    double p[3];
    knots->GetPoint(0, p);
    z = p[2];
    }

  if (knots->GetNumberOfPoints() < 2)
    {
    return;
    }

  double tmax, dmax;
  this->UpdateSpline(tmax, dmax);

  vtkIdType n = vtkMath::Floor(dmax) + 1;
  double delta = tmax/n;

  vtkSpline *xspline = this->SplineX;
  vtkSpline *yspline = this->SplineY;

  double t = tmax;
  if (xspline->GetClosed())
    {
    t = (n-1)*tmax/n;
    }
  else
    {
    n = n + 1;
    }

  t = 0;
  points->SetNumberOfPoints(0);
  for (vtkIdType i = 0; i < n; i++)
    {
    double p[3];

    p[0] = xspline->Evaluate(t);
    p[1] = yspline->Evaluate(t);
    p[2] = z;
    points->InsertNextPoint(p);

    t += delta;
    if (i == n-2)
      {
      t = 0;
      }
    }

  vtkIdType m = n;
  lines->Initialize();
  lines->InsertNextCell(m);
  for (vtkIdType i = 0; i < n; i++)
    {
    lines->InsertCellPoint(i);
    }
  if (m > n)
    {
    lines->InsertCellPoint(0);
    }

  this->ContourData->Modified();
}
