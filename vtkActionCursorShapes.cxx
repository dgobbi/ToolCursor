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

vtkCxxRevisionMacro(vtkActionCursorShapes, "$Revision: 1.1 $");
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

  data = this->MakePointerShape();
  this->AddShape("Pointer", data, VTK_SCURSOR_RGB);
  data->Delete();

  data = this->MakeCrosshairsShape();
  this->AddShape("Crosshairs", data, VTK_SCURSOR_RGB);
  data->Delete();

  data = this->MakeCrossShape(0);
  this->AddShape("Cross", data, 0);
  data->Delete();

  data = this->MakeCrossShape(1);
  this->AddShape("SplitCross", data, 0);
  data->Delete();

  data = this->MakeConeShape(0);
  this->AddShape("Cone", data, 0);
  data->Delete();

  data = this->MakeConeShape(1);
  this->AddShape("DualCone", data, 0);
  data->Delete();

  data = this->MakeSphereShape(0);
  this->AddShape("Sphere", data, 0);
  data->Delete();

  data = this->MakeSphereShape(1);
  this->AddShape("SplitSphere", data, 0);
  data->Delete();

  data = this->MakeMoverShape(0);
  this->AddShape("Mover", data, 0);
  data->Delete();

  data = this->MakeMoverShape(1);
  this->AddShape("Rocker", data, 0);
  data->Delete();

  data = this->MakePusherShape();
  this->AddShape("Pusher", data, VTK_SCURSOR_FLATX);
  data->Delete();

  data = this->MakeSpinnerShape();
  this->AddShape("Spinner", data, 0);
  data->Delete();
}

//----------------------------------------------------------------------------
vtkDataSet *vtkActionCursorShapes::MakePointerShape()
{
  vtkUnsignedCharArray *scalars = vtkUnsignedCharArray::New();
  scalars->SetNumberOfComponents(4);
  vtkPoints *points = vtkPoints::New();
  vtkCellArray *strips = vtkCellArray::New();  
  vtkCellArray *lines = vtkCellArray::New();  
  
  static unsigned char black[4] = {  0,   0,   0, 255};
  static unsigned char white[4] = {255, 255, 255, 255};

  static double hotspot[2] = { 0.5, -0.5 };
  static double coords[7][2] = {
    {  1,  -1 },
    {  1, -15 },
    {  4, -12 },
    {  8, -19 },
    { 10, -18 },
    {  7, -11 },
    { 11, -11 },
  };

  static vtkIdType stripIds[] = {
    3, 0, 1, 2,
    3, 0, 2, 5,
    3, 0, 5, 6,
    4, 2, 3, 5, 4,
  };

  static vtkIdType lineIds[] = {
    8, 7, 8, 9, 10, 11, 12, 13, 7,
    8, 14, 15, 16, 17, 18, 19, 20, 14,
  };

  // Add the points three times: white, then black, and black again
  for (int i = 0; i < 7; i++)
    {
    points->InsertNextPoint(coords[i][0] - hotspot[0],
                            coords[i][1] - hotspot[1], 0.0);
    scalars->InsertNextTupleValue(white);
    }

  for (int j = 0; j < 7; j++)
    {
    points->InsertNextPoint(coords[j][0] - hotspot[0],
                            coords[j][1] - hotspot[1], +0.1);
    scalars->InsertNextTupleValue(black);
    }

  for (int k = 0; k < 7; k++)
    {
    points->InsertNextPoint(coords[k][0] - hotspot[0],
                            coords[k][1] - hotspot[1], -0.1);
    scalars->InsertNextTupleValue(black);
    }

  // Make the strips
  strips->InsertNextCell(stripIds[0], &stripIds[1]);
  strips->InsertNextCell(stripIds[4], &stripIds[5]);
  strips->InsertNextCell(stripIds[8], &stripIds[9]);
  strips->InsertNextCell(stripIds[12], &stripIds[13]);

  // Make the lines
  lines->InsertNextCell(lineIds[0], &lineIds[1]);
  lines->InsertNextCell(lineIds[8], &lineIds[9]);

  vtkPolyData *data = vtkPolyData::New();
  data->SetPoints(points);
  points->Delete();
  data->SetStrips(strips);
  strips->Delete();
  data->SetLines(lines);
  lines->Delete();
  data->GetPointData()->SetScalars(scalars);
  scalars->Delete();

  return data;
}

//----------------------------------------------------------------------------
vtkDataSet *vtkActionCursorShapes::MakeCrosshairsShape()
{
  vtkUnsignedCharArray *scalars = vtkUnsignedCharArray::New();
  scalars->SetNumberOfComponents(4);
  vtkPoints *points = vtkPoints::New();
  vtkCellArray *strips = vtkCellArray::New();  
  vtkCellArray *lines = vtkCellArray::New();  
  
  static unsigned char black[4] = {  0,   0,   0, 255};
  static unsigned char white[4] = {255, 255, 255, 255};

  const double radius = 8;
  const double inner = 3;

  static double coords[8][2] = {
    { 0, -radius }, { 0, -inner }, { 0, +inner }, { 0, +radius },
    { -radius, 0 }, { -inner, 0 }, { +inner, 0 }, { +radius, 0 },
  };

  static double outCoords[16][2] = {
    { -1, -radius-1 }, { +1, -radius-1 }, { +1, -inner+1 }, { -1, -inner+1 }, 
    { +1, +radius+1 }, { -1, +radius+1 }, { -1, +inner-1 }, { +1, +inner-1 }, 
    { -radius-1, +1 }, { -radius-1, -1 }, { -inner+1, -1 }, { -inner+1, +1 }, 
    { +radius+1, -1 }, { +radius+1, +1 }, { +inner-1, +1 }, { +inner-1, -1 }, 
  };

  static vtkIdType toplineIds[] = {
    2, 0, 1,
    2, 2, 3,
    2, 4, 5,
    2, 6, 7,
  };

  static vtkIdType botlineIds[] = {
    2, 8, 9,
    2, 10, 11,
    2, 12, 13,
    2, 14, 15,
  };

  static vtkIdType outlineIds[] = {
    5, 16, 17, 18, 19, 16,
    5, 20, 21, 22, 23, 20,
    5, 24, 25, 26, 27, 24,
    5, 28, 29, 30, 31, 28,
  };

  static vtkIdType stripIds[] = {
    4, 16, 17, 19, 18,
    4, 20, 21, 23, 22,
    4, 24, 25, 27, 26,
    4, 28, 29, 31, 30,
  };

  for (int i = 0; i < 8; i++)
    {
    points->InsertNextPoint(coords[i][0]+0.5, coords[i][1]-0.5, +0.1);
    scalars->InsertNextTupleValue(black);
    }

  for (int j = 0; j < 8; j++)
    {
    points->InsertNextPoint(coords[j][0]+0.5, coords[j][1]-0.5, -0.1);
    scalars->InsertNextTupleValue(black);
    }

  for (int k = 0; k < 16; k++)
    {
    points->InsertNextPoint(outCoords[k][0]+0.5, outCoords[k][1]-0.5, 0);
    scalars->InsertNextTupleValue(white);
    }

  // Make the crosshairs
  lines->InsertNextCell(toplineIds[0], &toplineIds[1]);
  lines->InsertNextCell(toplineIds[3], &toplineIds[4]);
  lines->InsertNextCell(toplineIds[6], &toplineIds[7]);
  lines->InsertNextCell(toplineIds[9], &toplineIds[10]);

  // Make the crosshairs
  lines->InsertNextCell(botlineIds[0], &botlineIds[1]);
  lines->InsertNextCell(botlineIds[3], &botlineIds[4]);
  lines->InsertNextCell(botlineIds[6], &botlineIds[7]);
  lines->InsertNextCell(botlineIds[9], &botlineIds[10]);

  // Make the outline
  lines->InsertNextCell(outlineIds[0], &outlineIds[1]);
  lines->InsertNextCell(outlineIds[6], &outlineIds[7]);
  lines->InsertNextCell(outlineIds[12], &outlineIds[13]);
  lines->InsertNextCell(outlineIds[18], &outlineIds[19]);
  lines->InsertNextCell(outlineIds[24], &outlineIds[25]);

  // Fill the outline
  strips->InsertNextCell(stripIds[0], &stripIds[1]);
  strips->InsertNextCell(stripIds[5], &stripIds[6]);
  strips->InsertNextCell(stripIds[10], &stripIds[11]);
  strips->InsertNextCell(stripIds[15], &stripIds[16]);

  vtkPolyData *data = vtkPolyData::New();
  data->SetPoints(points);
  points->Delete();
  data->SetLines(lines);
  lines->Delete();
  data->SetStrips(strips);
  strips->Delete();
  data->GetPointData()->SetScalars(scalars);
  scalars->Delete();

  return data;
}

//----------------------------------------------------------------------------
vtkDataSet *vtkActionCursorShapes::MakeSphereShape(int dual)
{
  double pi = vtkMath::DoublePi();
  double radius = 5.0;
  int resolution = 9;

  vtkIdType *pointIds = new vtkIdType[4*(resolution+1)];

  vtkIntArray *scalars = vtkIntArray::New();
  vtkDoubleArray *normals = vtkDoubleArray::New();
  normals->SetNumberOfComponents(3);
  vtkPoints *points = vtkPoints::New();
  vtkCellArray *strips = vtkCellArray::New();
  vtkIdType nPoints = 0;
 
  int colorIndex = 0;

  for (int i = 0; i < 2; i++)
    {
    // The sign (i.e. for top or bottom) is stored in s
    double s = 1 - 2*colorIndex;

    // The unit position vector of the point is stored in v
    double v[3];
    v[0] = 0;
    v[1] = 0;
    v[2] = s;

    points->InsertNextPoint(radius*v[0], radius*v[1], radius*v[2]);
    normals->InsertNextTupleValue(v);
    scalars->InsertNextTupleValue(&colorIndex);
      
    int n = (resolution + 1)/2;
    int m = 2*resolution;

    for (int j = 1; j <= n; j++)
      {
      double phi = pi*j/resolution;
      double r = sin(phi);
      v[2] = cos(phi)*s;
      if (2*j >= resolution)
        {
        v[2] = 0;
        }

      for (int k = 0; k < m; k++)
        {
        double theta = pi*k/resolution;
        v[0] = r*cos(theta);
        v[1] = r*sin(theta)*s;
        points->InsertNextPoint(radius*v[0], radius*v[1], radius*v[2]);
        normals->InsertNextTupleValue(v);
        scalars->InsertNextTupleValue(&colorIndex);
        }
      }

    // Make the fan for the top
    pointIds[0] = nPoints++;
    for (int ii = 0; ii < (m-1); ii++)
      {
      pointIds[1] = nPoints + ii;
      pointIds[2] = nPoints + ii + 1;
      strips->InsertNextCell(3, pointIds);
      }
    pointIds[1] = nPoints + m - 1;
    pointIds[2] = nPoints;
    strips->InsertNextCell(3, pointIds);

    // Make the strips for the rest
    for (int jj = 1; jj < n; jj++)
      {
      for (int kk = 0; kk < m; kk++)
        {
        pointIds[2*kk] = nPoints + kk;
        pointIds[2*kk+1] = nPoints + kk + m;
        }
      pointIds[2*m] = nPoints;
      pointIds[2*m+1] = nPoints + m;
      strips->InsertNextCell(2*(m+1), pointIds);
      nPoints += m;
      }

    nPoints += m;

    if (dual)
      {
      colorIndex = 1;
      }
    }

  delete [] pointIds;

  vtkPolyData *data = vtkPolyData::New();
  data->SetPoints(points);
  points->Delete();
  data->SetStrips(strips);
  strips->Delete();
  data->GetPointData()->SetScalars(scalars);
  scalars->Delete();
  data->GetPointData()->SetNormals(normals);
  normals->Delete();

  return data;
}

//----------------------------------------------------------------------------
vtkDataSet *vtkActionCursorShapes::MakeConeShape(int dual)
{
  double pi = vtkMath::DoublePi();
  double radius = 8.0;
  double height = 15.0;
  int resolution = 20;

  vtkIdType *pointIds = new vtkIdType[2*(resolution+1)];

  vtkIntArray *scalars = vtkIntArray::New();
  vtkDoubleArray *normals = vtkDoubleArray::New();
  normals->SetNumberOfComponents(3);
  vtkPoints *points = vtkPoints::New();
  vtkCellArray *strips = vtkCellArray::New();
  vtkIdType nPoints = 0;
 
  int sides = (dual ? 2 : 1);

  for (int colorIndex = 0; colorIndex < sides; colorIndex++)
    {
    // The sign (i.e. for top or bottom) is stored in s
    double s = 1 - 2*colorIndex;

    // The length of the side of the cone
    double l = sqrt(radius*radius + height*height);
    double f1 = radius/l;
    double f2 = height/l;

    // The unit normal vector
    double v[3];

    // The point of the cone
    for (int i = 0; i < 2; i++)
      {
      double r = radius*i;
      double z = height*i;
      double offset = 0.5*(1 - i);

      for (int j = 0; j < resolution; j++)
        {
        double theta = 2*pi*(j + offset)/resolution;
        double ct = cos(theta);
        double st = sin(theta);
        v[0] = f2*ct;
        v[1] = f2*st*s;
        v[2] = -f1*s;
        points->InsertNextPoint(r*ct, r*st*s, z*s);
        normals->InsertNextTupleValue(v);
        scalars->InsertNextTupleValue(&colorIndex);
        }
      }

    // The base of the cone
    v[0] = 0;
    v[1] = 0;
    v[2] = s;
    points->InsertNextPoint(0, 0, height*s);
    normals->InsertNextTupleValue(v);
    scalars->InsertNextTupleValue(&colorIndex); 

    for (int k = 0; k < resolution; k++)
      {
      double theta = 2*pi*k/resolution;
      points->InsertNextPoint(radius*cos(theta), radius*sin(theta)*s, height*s);
      normals->InsertNextTupleValue(v);
      scalars->InsertNextTupleValue(&colorIndex);
      }

    // Make the fan for the top
    for (int ii = 0; ii < (resolution-1); ii++)
      {
      pointIds[0] = nPoints + ii;
      pointIds[1] = nPoints + ii + resolution + 1;
      pointIds[2] = nPoints + ii + resolution;
      strips->InsertNextCell(3, pointIds);
      }
    pointIds[0] = nPoints + 2*resolution - 1;
    pointIds[1] = nPoints + resolution - 1;
    pointIds[2] = nPoints + resolution;
    strips->InsertNextCell(3, pointIds);
    nPoints += 2*resolution;
 
    // Make the fan for the base
    pointIds[0] = nPoints++;
    for (int jj = 0; jj < (resolution-1); jj++)
      {
      pointIds[1] = nPoints + jj;
      pointIds[2] = nPoints + jj + 1;
      strips->InsertNextCell(3, pointIds);
      }
    pointIds[1] = nPoints + resolution - 1;
    pointIds[2] = nPoints;
    strips->InsertNextCell(3, pointIds);
    nPoints += resolution;
    }

  delete [] pointIds;

  vtkPolyData *data = vtkPolyData::New();
  data->SetPoints(points);
  points->Delete();
  data->SetStrips(strips);
  strips->Delete();
  data->GetPointData()->SetScalars(scalars);
  scalars->Delete();
  data->GetPointData()->SetNormals(normals);
  normals->Delete();

  return data;
}

//----------------------------------------------------------------------------
vtkPolyData *vtkActionCursorShapes::MakeWarpedArrow(double warpX, double warpY,
                                               double warpZ, double warpScale)
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

  vtkWarpTo *warp = vtkWarpTo::New();
  warp->SetInputConnection(reverse->GetOutputPort());
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
  reverse->Delete();
  contour->Delete();
  modeller->Delete();

  return data;
}

//----------------------------------------------------------------------------
vtkDataSet *vtkActionCursorShapes::MakeMoverShape(int warped)
{
  vtkPolyData *leafData = vtkActionCursorShapes::MakeWarpedArrow(
    10, 0, -10, (warped ? 1.0 : 0.5));
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
vtkDataSet *vtkActionCursorShapes::MakePusherShape()
{
  vtkPolyData *leafData = vtkActionCursorShapes::MakeWarpedArrow(
    10.0, 0.0, 10.0, 0.0);
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
vtkDataSet *vtkActionCursorShapes::MakeSpinnerShape()
{
  vtkPolyData *leafData = vtkActionCursorShapes::MakeWarpedArrow(
    10.0, 0.0, -5.0, 1.0);
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

//----------------------------------------------------------------------------
vtkDataSet *vtkActionCursorShapes::MakeCrossShape(int dual)
{
  double radius = 10.0;
  double inner = 3.5;
  double thickness = 2.0;

  double xmin = inner;
  double xmax = radius;
  double ymin = -thickness*0.5;
  double ymax = +thickness*0.5;
  double zmin = 0;
  double zmax = thickness*0.5;

  vtkIntArray *scalars = vtkIntArray::New();
  vtkPoints *points = vtkPoints::New();
  vtkCellArray *strips = vtkCellArray::New();
  vtkIdType nPoints = 0;
  
  int colorIndex = 0;

  for (int i = 0; i < 2; i++)
    {
    for (int j = 0; j < 4; j++)
      {
      double z = zmax;
      for (int k = 0; k < 2; k++)
        {
        points->InsertNextPoint(xmin, ymin, z);
        points->InsertNextPoint(xmin, ymax, z);
        points->InsertNextPoint(xmax, ymax, z);
        points->InsertNextPoint(xmax, ymin, z);
        scalars->InsertNextTupleValue(&colorIndex);
        scalars->InsertNextTupleValue(&colorIndex);
        scalars->InsertNextTupleValue(&colorIndex);
        scalars->InsertNextTupleValue(&colorIndex);
        z = zmin;
        }
   
      static vtkIdType rectIds[5][4] = { { 1, 0, 2, 3 },
                                         { 4, 0, 5, 1 },
                                         { 5, 1, 6, 2 },
                                         { 6, 2, 7, 3 },
                                         { 7, 3, 4, 0 } };
      vtkIdType pointIds[4];
      for (int ii = 0; ii < 5; ii++)
        {
        for (int jj = 0; jj < 4; jj++)
          {
          pointIds[jj] = rectIds[ii][jj]+nPoints;
          }
        strips->InsertNextCell(4, pointIds);
        }
      nPoints += 8;

      // do a rotation of 90 degrees for next piece
      double tmp1 = ymin;
      double tmp2 = ymax;
      ymin = -xmax;
      ymax = -xmin;
      xmin = tmp1;
      xmax = tmp2;
      }

    // do the other side
    zmin = -zmin;
    zmax = -zmax;
    xmin = -xmin;
    xmax = -xmax;

    if (dual)
      {
      colorIndex = 1;
      }
    }

  vtkPolyData *data = vtkPolyData::New();
  data->SetPoints(points);
  points->Delete();
  data->SetStrips(strips);
  strips->Delete();
  data->GetPointData()->SetScalars(scalars);
  scalars->Delete();

  return data;
}

