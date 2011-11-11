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
#include "vtkGlyph3D.h"
#include "vtkCellLocator.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
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

  this->ROIData = vtkROIContourData::New();

  this->ROIDataToPointSet = vtkROIContourDataToPolyDataFilter::New();
  this->ROIDataToPointSet->SetInput(this->ROIData);
  this->ROIDataToPolyData = vtkROIContourDataToPolyDataFilter::New();
  this->ROIDataToPolyData->SetInput(this->ROIData);
  this->ROIDataToPolyData->SubdivisionOn();

  this->Glyph3D = vtkGlyph3D::New();
  this->Glyph3D->SetColorModeToColorByScalar();
  this->Glyph3D->SetScaleModeToDataScalingOff();
  this->Glyph3D->SetScaleFactor(0.2);
  this->Glyph3D->SetInputConnection(this->ROIDataToPointSet->GetOutputPort());
  this->Glyph3D->SetSource(vtkPolyData::SafeDownCast(
    shapes->GetShapeData("Sphere")));

  this->GlyphMapper = vtkDataSetMapper::New();
  this->GlyphMapper->SetInputConnection(this->Glyph3D->GetOutputPort());

  this->GlyphActor = vtkActor::New();
  this->GlyphActor->PickableOff();
  this->GlyphActor->SetMapper(this->GlyphMapper);
  this->GlyphActor->GetProperty()->SetColor(1,0,0);
  this->GlyphActor->GetProperty()->LightingOff();
  this->GlyphActor->SetUserMatrix(this->Matrix);

  this->ContourMapper = vtkDataSetMapper::New();
  this->ContourMapper->SetInputConnection(this->ROIDataToPolyData->GetOutputPort());

  this->ContourActor = vtkActor::New();
  this->ContourActor->PickableOff();
  this->ContourActor->SetMapper(this->ContourMapper);
  this->ContourActor->GetProperty()->SetColor(1,0,0);
  this->ContourActor->GetProperty()->LightingOff();
  //this->ContourActor->GetProperty()->SetLineStipplePattern(0xe0e0e0e0);
  this->ContourActor->GetProperty()->SetLineStipplePattern(0xfcfcfcfc);
  this->ContourActor->SetUserMatrix(this->Matrix);

  this->CellLocator = vtkCellLocator::New();

  shapes->Delete();

  this->CurrentPointId = -1;
  this->CurrentContourId = -1;
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

  this->CellLocator->Delete();
  this->Matrix->Delete();
  this->ROIDataToPointSet->Delete();
  this->ROIDataToPolyData->Delete();
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

  // Tolerance for point selection
  double tol = 3;

  // Get the spline data and create a locator
  this->ROIDataToPolyData->Update();
  vtkPolyData *contourData =
    vtkPolyData::SafeDownCast(this->ROIDataToPolyData->GetOutput());
  vtkCellLocator *cellLocator = 0;
  if (contourData->GetNumberOfCells() > 0)
    {
    cellLocator = this->CellLocator;
    cellLocator->SetDataSet(contourData);
    cellLocator->BuildLocator();
    }

  // Get the contour data
  vtkROIContourData *data = this->ROIData;
  int numContours = data->GetNumberOfContours();
  this->CurrentPointId = -1;
  this->CurrentContourId = -1;
  double tol2 = tol*tol;
  for (int ic = 0; ic < numContours; ic++)
    {
    // Check if mouse is over a point
    vtkPoints *points = data->GetContourPoints(ic);
    vtkIdType n = points->GetNumberOfPoints();
    for (vtkIdType i = 0; i < n; i++)
      {
      double p[3];
      points->GetPoint(i, p);
      double d2 = vtkMath::Distance2BetweenPoints(position, p);
      if (d2 <= tol2)
        {
        tol2 = d2;
        this->CurrentPointId = i;
        this->CurrentContourId = ic;
        this->InitialPointPosition[0] = p[0];
        this->InitialPointPosition[1] = p[1];
        this->InitialPointPosition[2] = p[2];
        break;
        }
      }
    }

  if (this->CurrentPointId == 0)
    {
    // Clicked on first point: close the contour
    this->ROIData->SetContourType(
      this->CurrentContourId, vtkROIContourData::CLOSED_PLANAR);
    }
  else if (this->CurrentPointId < 0)
    {
    vtkIdType cellId;
    int subId;
    double d2;
    double p[3];
    if (cellLocator &&
        (cellLocator->FindClosestPointWithinRadius(
         position, tol, p, cellId, subId, d2)))
      {
      vtkIntArray *contourIds = vtkIntArray::SafeDownCast(
        contourData->GetCellData()->GetArray("Labels"));
      vtkIntArray *contourSubIds = vtkIntArray::SafeDownCast(
        contourData->GetPointData()->GetArray("SubIds"));

      // Get the cell for this contour
      vtkCellArray *lines = contourData->GetLines();
      vtkCellArray *verts = contourData->GetVerts();
      if (verts)
        {
        // verts come before lines, so adjust cell Id
        cellId -= verts->GetNumberOfCells();
        }
      vtkIdType npts, *pts;
      vtkIdType l = 0;
      vtkIdType icell = cellId + 1;
      do
        {
        lines->GetCell(l, npts, pts);
        l += npts + 1;
        }
      while (--icell > 0);

      // Get the pointId for insertion, and get the contourId
      subId = contourSubIds->GetValue(pts[subId + 1]);
      int contourId = contourIds->GetValue(cellId);
      vtkPoints *points = data->GetContourPoints(contourId);

      // Insert the point at the correct position
      int m = static_cast<int>(points->InsertNextPoint(p));
      for (int i = m; i > subId + 1; --i)
        {
        double ptmp[3];
        points->GetPoint(i-1, ptmp);
        points->SetPoint(i, ptmp); 
        }
      points->SetPoint(subId + 1, p);

      this->CurrentContourId = contourId;
      this->CurrentPointId = subId + 1;
      this->InitialPointPosition[0] = p[0];
      this->InitialPointPosition[1] = p[1];
      this->InitialPointPosition[2] = p[2];
      }

    if (this->CurrentContourId < 0)
      {
      // Check if there is an open contour
      for (int i = 0; i < numContours; i++)
        {
        if (data->GetContourType(i) == vtkROIContourData::OPEN_PLANAR)
          {
          this->CurrentContourId = i;
          }
        }
      }

    if (this->CurrentContourId < 0)
      {
      // Create a new contour
      vtkPoints *points = vtkPoints::New(VTK_DOUBLE);
      points->InsertNextPoint(position);
      data->SetNumberOfContours(numContours+1);
      data->SetContourPoints(numContours, points);
      data->SetContourType(numContours, vtkROIContourData::OPEN_PLANAR);
      }
    else if (this->CurrentPointId < 0)
      {
      // Add a new point to the currently unclosed contour
      vtkPoints *points = data->GetContourPoints(this->CurrentContourId);
      this->CurrentPointId = points->InsertNextPoint(position);
      this->InitialPointPosition[0] = position[0];
      this->InitialPointPosition[1] = position[1];
      this->InitialPointPosition[2] = position[2];
      }

    this->ROIData->Modified();
    cursor->GetRenderer()->AddViewProp(this->GlyphActor);
    cursor->GetRenderer()->AddViewProp(this->ContourActor);
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
  double dz = position[2] - p0[2];

  vtkPoints *points = 0;
  if (this->CurrentContourId >= 0)
    {
    points = this->ROIData->GetContourPoints(this->CurrentContourId);
    }

  if (points && this->CurrentPointId >= 0)
    {
    double p[3];
    p[0] = this->InitialPointPosition[0] + dx;
    p[1] = this->InitialPointPosition[1] + dy;
    p[2] = this->InitialPointPosition[2] + dz;
    points->SetPoint(this->CurrentPointId, p);

    this->ROIData->Modified();
    }
}
