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

  this->ROIData = vtkROIContourData::New();

  this->ROIDataToPointSet = vtkROIContourDataToPolyDataFilter::New();
  this->ROIDataToPointSet->SetInput(this->ROIData);
  this->ROIDataToPolyData = vtkROIContourDataToPolyDataFilter::New();
  this->ROIDataToPolyData->SetInput(this->ROIData);
  this->ROIDataToPolyData->SubdivisionOn();

  this->Glyph3D = vtkGlyph3D::New();
  this->Glyph3D->SetColorModeToColorByScalar();
  this->Glyph3D->SetScaleFactor(0.7);
  this->Glyph3D->SetInputConnection(this->ROIDataToPointSet->GetOutputPort());
  this->Glyph3D->SetSource(vtkPolyData::SafeDownCast(
    shapes->GetShapeData("Sphere")));

  this->GlyphMapper = vtkDataSetMapper::New();
  this->GlyphMapper->SetInputConnection(this->Glyph3D->GetOutputPort());

  this->GlyphActor = vtkActor::New();
  this->GlyphActor->PickableOff();
  this->GlyphActor->SetMapper(this->GlyphMapper);
  this->GlyphActor->GetProperty()->SetColor(1,0,0);
  this->GlyphActor->SetUserMatrix(this->Matrix);

  this->ContourMapper = vtkDataSetMapper::New();
  this->ContourMapper->SetInputConnection(this->ROIDataToPolyData->GetOutputPort());

  this->ContourActor = vtkActor::New();
  this->ContourActor->PickableOff();
  this->ContourActor->SetMapper(this->ContourMapper);
  this->ContourActor->GetProperty()->SetColor(1,0,0);
  //this->ContourActor->GetProperty()->SetLineStipplePattern(0xe0e0e0e0);
  this->ContourActor->GetProperty()->SetLineStipplePattern(0xfcfcfcfc);
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
    }
}
