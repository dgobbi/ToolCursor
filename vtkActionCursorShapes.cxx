/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkActionCursorShapes.cxx,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkActionCursorShapes.h"
#include "vtkObjectFactory.h"

#include "vtkSurfaceCursor.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkDataArray.h"
#include "vtkDoubleArray.h"
#include "vtkIntArray.h"
#include "vtkPointData.h"
#include "vtkMath.h"

#include "vtkImplicitModeller.h"
#include "vtkContourFilter.h"
#include "vtkStripper.h"
#include "vtkPolyDataNormals.h"
#include "vtkReverseSense.h"
#include "vtkWarpTo.h"
#include "vtkTransform.h"

vtkCxxRevisionMacro(vtkActionCursorShapes, "$Revision: 1.3 $");
vtkStandardNewMacro(vtkActionCursorShapes);

//----------------------------------------------------------------------------
vtkActionCursorShapes::vtkActionCursorShapes()
{
  this->MakeShapes();
}

//----------------------------------------------------------------------------
vtkActionCursorShapes::~vtkActionCursorShapes()
{
}

//----------------------------------------------------------------------------
void vtkActionCursorShapes::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

//----------------------------------------------------------------------------
void vtkActionCursorShapes::MakeShapes()
{
  vtkDataSet *data;

  vtkPolyData *arrow = this->MakeArrow();

  data = this->MakeMoverShape(arrow, 0);
  this->AddShape("Mover", data, 0);
  data->Delete();

  data = this->MakeMoverShape(arrow, 1);
  this->AddShape("Rotator", data, 0);
  data->Delete();

  data = this->MakePusherShape(arrow);
  this->AddShape("Pusher", data, VTK_SCURSOR_FLATX);
  data->Delete();

  data = this->MakeSpinnerShape(arrow);
  this->AddShape("Spinner", data, 0);
  data->Delete();

  arrow->Delete();
}

//----------------------------------------------------------------------------
vtkPolyData *vtkActionCursorShapes::MakeArrow()
{
  vtkPoints *points = vtkPoints::New();
  vtkCellArray *polys = vtkCellArray::New();

  static double coords[7][3] = {
    { 4, 2, 0 }, { 15, 2, 0 }, { 14, 8, 0 }, { 24, 0.01, 0 },
    { 14, -8, 0 }, { 15, -2, 0 }, { 4, -2, 0 },
  };

  static vtkIdType polyIds[] = {
    7, 0, 1, 2, 3, 4, 5, 6,
  };

  for (int i = 0; i < 7; i++)
    {
    points->InsertNextPoint(coords[i]);
    }
  polys->InsertNextCell(polyIds[0], &polyIds[1]);

  vtkPolyData *arrow = vtkPolyData::New();
  arrow->SetPoints(points);
  points->Delete();
  arrow->SetPolys(polys);
  polys->Delete();

  vtkImplicitModeller *modeller = vtkImplicitModeller::New();
  modeller->SetInput(arrow);
  modeller->SetSampleDimensions(32, 16, 8);

  vtkContourFilter *contour = vtkContourFilter::New();
  contour->SetInputConnection(modeller->GetOutputPort());
  contour->SetValue(0, 0.5);

  // The image is inside-out, and so is the contour
  vtkReverseSense *reverse = vtkReverseSense::New();
  reverse->SetInputConnection(contour->GetOutputPort());
  reverse->ReverseCellsOn();

  reverse->Update();

  vtkPolyData *data = vtkPolyData::New();
  data->DeepCopy(reverse->GetOutput());

  reverse->Delete();
  contour->Delete();
  modeller->Delete();

  return data;
}

//----------------------------------------------------------------------------
vtkPolyData *vtkActionCursorShapes::MakeWarpedArrow(vtkPolyData *arrow,
  double warpX, double warpY, double warpZ, double warpScale)
{
  vtkWarpTo *warp = vtkWarpTo::New();
  warp->SetInput(arrow);
  warp->SetPosition(warpX, warpY, warpZ);
  warp->SetScaleFactor(warpScale);
  warp->AbsoluteOn();

  vtkPolyDataNormals *polyNormals = vtkPolyDataNormals::New();
  polyNormals->SetInputConnection(warp->GetOutputPort());
  polyNormals->SplittingOn();

  vtkStripper *stripper = vtkStripper::New();
  stripper->SetInputConnection(polyNormals->GetOutputPort());

  stripper->Update();

  vtkPolyData *data = vtkPolyData::New();
  data->DeepCopy(stripper->GetOutput());

  polyNormals->Delete();
  stripper->Delete();
  warp->Delete();

  return data;
}

//----------------------------------------------------------------------------
vtkDataSet *vtkActionCursorShapes::MakeMoverShape(
  vtkPolyData *arrow, int warped)
{
  vtkPolyData *leafData = vtkActionCursorShapes::MakeWarpedArrow(
    arrow, 10, 0, -10, (warped ? 1.0 : 0.5));
  vtkPoints *leafPoints = leafData->GetPoints();
  vtkCellArray *leafStrips = leafData->GetStrips();
  vtkDataArray *leafNormals = leafData->GetPointData()->GetNormals();

  vtkPolyData *data = vtkPolyData::New();
  vtkPoints *points = vtkPoints::New();
  vtkDoubleArray *normals = vtkDoubleArray::New();
  normals->SetNumberOfComponents(3);
  vtkCellArray *strips = vtkCellArray::New();
  vtkIntArray *scalars = vtkIntArray::New();

  vtkTransform *transform = vtkTransform::New();
  transform->PostMultiply();

  static double rotate90[16] = {
    0,-1, 0, 0,
    1, 0, 0, 0,
    0, 0, 1, 0,
    0, 0, 0, 1,
  };

  if (warped)
    {
    transform->RotateY(36);
    transform->Scale(1.5, 1.0, 1.0);
    transform->Translate(0, 0, 14);
    }

  int color = 0;

  for (int j = 0; j < 4; j++)
    {
    vtkIdType nPoints = points->GetNumberOfPoints();
    transform->TransformPoints(leafPoints, points);
    transform->TransformNormals(leafNormals, normals);
    vtkIdType npts;
    vtkIdType *pts;
    leafStrips->InitTraversal();
    while (leafStrips->GetNextCell(npts, pts))
      {
      strips->InsertNextCell(npts);
      for (vtkIdType k = 0; k < npts; k++)
        {
        strips->InsertCellPoint(pts[k] + nPoints);
        }
      }
    vtkIdType n = leafPoints->GetNumberOfPoints();
    for (int ii = 0; ii < n; ii++)
      {
      scalars->InsertNextTupleValue(&color);
      }
    transform->Concatenate(rotate90);
    color = !color;
    }
  transform->Delete();

  leafData->Delete();

  data->SetPoints(points);
  points->Delete();
  data->GetPointData()->SetNormals(normals);
  normals->Delete();
  data->GetPointData()->SetScalars(scalars);
  scalars->Delete();
  data->SetStrips(strips);
  strips->Delete();

  return data;
}

//----------------------------------------------------------------------------
vtkDataSet *vtkActionCursorShapes::MakePusherShape(vtkPolyData *arrow)
{
  vtkPolyData *leafData = vtkActionCursorShapes::MakeWarpedArrow(
    arrow, 10.0, 0.0, 10.0, 0.0);
  vtkPoints *leafPoints = leafData->GetPoints();
  vtkCellArray *leafStrips = leafData->GetStrips();
  vtkDataArray *leafNormals = leafData->GetPointData()->GetNormals();

  vtkPolyData *data = vtkPolyData::New();
  vtkPoints *points = vtkPoints::New();
  vtkDoubleArray *normals = vtkDoubleArray::New();
  normals->SetNumberOfComponents(3);
  vtkCellArray *strips = vtkCellArray::New();
  vtkCellArray *lines = vtkCellArray::New();
  vtkIntArray *scalars = vtkIntArray::New();

  vtkTransform *transform = vtkTransform::New();
  transform->PostMultiply();

  static double rotate90[16] = {
    0, 0, 1, 0,
    0, 1, 0, 0,
   -1, 0, 0, 0,
    0, 0, 0, 1,
  };

  transform->Concatenate(rotate90);
  transform->Scale(1.0, 1.0, 1.0);
  transform->Translate(0, 0, 24);

  int color = 0;

  for (int i = 0; i < 2; i++)
    {
    vtkIdType nPoints = points->GetNumberOfPoints();
    transform->TransformPoints(leafPoints, points);
    transform->TransformNormals(leafNormals, normals);
    vtkIdType npts;
    vtkIdType *pts;
    leafStrips->InitTraversal();
    while (leafStrips->GetNextCell(npts, pts))
      {
      strips->InsertNextCell(npts);
      for (vtkIdType j = 0; j < npts; j++)
        {
        strips->InsertCellPoint(pts[j] + nPoints);
        }
      }
    vtkIdType nn = leafPoints->GetNumberOfPoints();
    for (int ii = 0; ii < nn; ii++)
      {
      scalars->InsertNextTupleValue(&color);
      }
    // This transform puts the arrows tail-to-tail, instead of tip-to-tip 
    transform->Translate(0, 0, -44);
    transform->Scale(-1,1,-1);
    color = !color;
    }
  transform->Delete();

  leafData->Delete();

  // Make a ring for when arrow viewed tail-on
  const double lineRadius = 8;
  const int lineResolution = 24; // must be divisible by 8
  int polylen = lineResolution/8 + 1;
  double normal[3];
  normal[0] = normal[1] = 0.0;
  normal[2] = 1.0;
  for (int j = 0; j < 2; j++)
    {
    color = 0;
    for (int k = 0; k < 8; k++)
      {
      vtkIdType nPoints = points->GetNumberOfPoints();
      lines->InsertNextCell(polylen);
      for (int ii = 0; ii < polylen; ii++)
        {
        double angle = 2*vtkMath::DoublePi()*(k*(polylen-1)+ii)/lineResolution;
        points->InsertNextPoint(lineRadius*cos(angle),
                                lineRadius*sin(angle),
                                0.1*(1 - 2*j));
        scalars->InsertNextTupleValue(&color);
        normals->InsertNextTupleValue(normal);
        lines->InsertCellPoint(ii + nPoints);
        }
      normal[2] = -normal[2];
      color = !color;
      }
    }

  data->SetPoints(points);
  points->Delete();
  data->GetPointData()->SetNormals(normals);
  normals->Delete();
  data->GetPointData()->SetScalars(scalars);
  scalars->Delete();
  data->SetStrips(strips);
  strips->Delete();
  data->SetLines(lines);
  lines->Delete();

  return data;
}

//----------------------------------------------------------------------------
vtkDataSet *vtkActionCursorShapes::MakeSpinnerShape(vtkPolyData *arrow)
{
  vtkPolyData *leafData = vtkActionCursorShapes::MakeWarpedArrow(
    arrow, 10.0, 0.0, -5.0, 1.0);
  vtkPoints *leafPoints = leafData->GetPoints();
  vtkCellArray *leafStrips = leafData->GetStrips();
  vtkDataArray *leafNormals = leafData->GetPointData()->GetNormals();

  vtkPolyData *data = vtkPolyData::New();
  vtkPoints *points = vtkPoints::New();
  vtkDoubleArray *normals = vtkDoubleArray::New();
  normals->SetNumberOfComponents(3);
  vtkCellArray *strips = vtkCellArray::New();
  vtkIntArray *scalars = vtkIntArray::New();

  vtkTransform *transform = vtkTransform::New();
  transform->PostMultiply();

  static double rotate90[16] = {
    1, 0, 0, 0,
    0, 0, 1, 0,
    0,-1, 0, 0,
    0, 0, 0, 1,
  };

  transform->Concatenate(rotate90);
  transform->Scale(3.0, 3.0, 1.5);
  transform->Translate(-30, 16, 3.75);

  int color = 0;

  for (int j = 0; j < 2; j++)
    {
    vtkIdType nPoints = points->GetNumberOfPoints();
    transform->TransformPoints(leafPoints, points);
    transform->TransformNormals(leafNormals, normals);
    vtkIdType npts;
    vtkIdType *pts;
    leafStrips->InitTraversal();
    while (leafStrips->GetNextCell(npts, pts))
      {
      strips->InsertNextCell(npts);
      for (vtkIdType k = 0; k < npts; k++)
        {
        strips->InsertCellPoint(pts[k] + nPoints);
        }
      }
    vtkIdType n = leafPoints->GetNumberOfPoints();
    for (int ii = 0; ii < n; ii++)
      {
      scalars->InsertNextTupleValue(&color);
      }
    transform->Scale(-1, -1, 1);
    color = !color;
    }
  transform->Delete();

  leafData->Delete();

  data->SetPoints(points);
  points->Delete();
  data->GetPointData()->SetNormals(normals);
  normals->Delete();
  data->GetPointData()->SetScalars(scalars);
  scalars->Delete();
  data->SetStrips(strips);
  strips->Delete();

  return data;
}
