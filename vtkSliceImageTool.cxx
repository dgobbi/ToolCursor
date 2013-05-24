/*=========================================================================

  Program:   ToolCursor
  Module:    vtkSliceImageTool.cxx

  Copyright (c) 2010 David Gobbi
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkSliceImageTool.h"
#include "vtkObjectFactory.h"

#include "vtkToolCursor.h"
#include "vtkCamera.h"
#include "vtkRenderer.h"
#include "vtkImageMapper3D.h"
#include "vtkImageData.h"
#include "vtkMatrix4x4.h"
#include "vtkMath.h"

#include "vtkVolumePicker.h"

vtkStandardNewMacro(vtkSliceImageTool);

//----------------------------------------------------------------------------
vtkSliceImageTool::vtkSliceImageTool()
{
  this->JumpToNearestSlice = 0;
}

//----------------------------------------------------------------------------
vtkSliceImageTool::~vtkSliceImageTool()
{
}

//----------------------------------------------------------------------------
void vtkSliceImageTool::SetJumpToNearestSlice(int val)
{
  if (val != this->JumpToNearestSlice)
    {
    this->JumpToNearestSlice = val;
    this->Modified();
    }
}

//----------------------------------------------------------------------------
void vtkSliceImageTool::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

//----------------------------------------------------------------------------
void vtkSliceImageTool::StartAction()
{
  this->Superclass::StartAction();

  vtkToolCursor *cursor = this->GetToolCursor();
  vtkCamera *camera = cursor->GetRenderer()->GetActiveCamera();

  this->StartDistance = camera->GetDistance();

    // code for handling the mouse wheel interaction
  if ((cursor->GetModifier() & VTK_TOOL_WHEEL_MASK) != 0)
    {
    int delta = 0;
    if ((cursor->GetModifier() & VTK_TOOL_WHEEL_BWD) != 0)
      {
      delta = -1;
      }
    else if ((cursor->GetModifier() & VTK_TOOL_WHEEL_FWD) != 0)
      {
      delta = 1;
      }
    this->AdvanceSlice(delta);
    }
}

//----------------------------------------------------------------------------
void vtkSliceImageTool::StopAction()
{
  this->Superclass::StopAction();
}

//----------------------------------------------------------------------------
void vtkSliceImageTool::DoAction()
{
  this->Superclass::DoAction();

  vtkToolCursor *cursor = this->GetToolCursor();
  vtkCamera *camera = cursor->GetRenderer()->GetActiveCamera();

  // Get the display position.
  double x0, y0, x, y;
  this->GetStartDisplayPosition(x0, y0);
  this->GetDisplayPosition(x, y);

  // Get viewport height at the current depth
  double height = 1;
  if (camera->GetParallelProjection())
    {
    height = camera->GetParallelScale();
    }
  else
    {
    double angle = vtkMath::RadiansFromDegrees(camera->GetViewAngle());
    height = 2*this->StartDistance*sin(angle*0.5);
    }

  // Get the viewport size in pixels
  vtkRenderer *renderer = cursor->GetRenderer();
  int *size = renderer->GetSize();

  // Get the vertical offset
  double dy = y - y0;
  double distance = this->StartDistance + dy*height/size[1];

  double focalPoint[3], position[3];
  camera->GetFocalPoint(focalPoint);
  camera->GetPosition(position);

  double vector[3];
  vector[0] = focalPoint[0] - position[0];
  vector[1] = focalPoint[1] - position[1];
  vector[2] = focalPoint[2] - position[2];

  vtkMath::Normalize(vector);

  focalPoint[0] = position[0] + distance*vector[0];
  focalPoint[1] = position[1] + distance*vector[1];
  focalPoint[2] = position[2] + distance*vector[2];

  vtkImageData *data = 0;

  if (this->JumpToNearestSlice &&
      this->CurrentImageMatrix &&
      this->CurrentImageMapper &&
      (data = vtkImageData::SafeDownCast(this->CurrentImageMapper->GetInput())))
    {
    double point[4];
    point[0] = focalPoint[0];
    point[1] = focalPoint[1];
    point[2] = focalPoint[2];
    point[3] = 1.0;

    double normal[4];
    normal[0] = vector[0];
    normal[1] = vector[1];
    normal[2] = vector[2];
    normal[3] = -vtkMath::Dot(point, normal);

    // convert normal to data coordinates
    double worldToData[16];
    vtkMatrix4x4 *dataToWorld = this->CurrentImageMatrix;
    vtkMatrix4x4::Transpose(*dataToWorld->Element, worldToData);
    vtkMatrix4x4::MultiplyPoint(worldToData, normal, normal);

    // find the slice orientation from the normal
    int k = 0;
    double maxsq = 0;
    double sumsq = 0;
    for (int i = 0; i < 3; i++)
      {
      double tmpsq = normal[i]*normal[i];
      sumsq += tmpsq;
      if (tmpsq > maxsq)
        {
        maxsq = tmpsq;
        k = i;
        }
      }

    // if the slice is not oblique
    if ((1.0 - maxsq/sumsq) < 1e-12)
      {
      // get the point in data coordinates
      vtkMatrix4x4::Invert(*dataToWorld->Element, worldToData);
      vtkMatrix4x4::MultiplyPoint(worldToData, point, point);

      double *spacing = data->GetSpacing();
      double *origin = data->GetOrigin();

      // set the point to lie exactly on a slice
      double z = (point[k] - origin[k])/spacing[k];
      if (z > VTK_INT_MIN && z < VTK_INT_MAX)
        {
        int j = vtkMath::Floor(z + 0.5);
        point[k] = j*spacing[k] + origin[k];
        }

      // convert back to world coordinates
      dataToWorld->MultiplyPoint(point, point);

      if (point[3] != 0)
        {
        focalPoint[0] = point[0]/point[3];
        focalPoint[1] = point[1]/point[3];
        focalPoint[2] = point[2]/point[3];
        }
      }
    }

  camera->SetFocalPoint(focalPoint);
}

//----------------------------------------------------------------------------
void vtkSliceImageTool::AdvanceSlice(int delta)
{
  vtkToolCursor *cursor = this->GetToolCursor();
  vtkCamera *camera = cursor->GetRenderer()->GetActiveCamera();

  // Get the camera information
  double point[4], normal[4], position[3];
  camera->GetFocalPoint(point);
  camera->GetPosition(position);
  normal[0] = point[0] - position[0];
  normal[1] = point[1] - position[1];
  normal[2] = point[2] - position[2];
  vtkMath::Normalize(normal);
  normal[3] = -vtkMath::Dot(point, normal);

  // Convert the normal to data coordinates
  if (this->CurrentImageMatrix)
    {
    double matrix[16];
    vtkMatrix4x4::Transpose(*this->CurrentImageMatrix->Element, matrix);
    vtkMatrix4x4::MultiplyPoint(matrix, normal, normal);
    }

  // Get the image information
  vtkImageData *data = 0;
  if (this->CurrentImageMapper)
    {
    data = vtkImageData::SafeDownCast(this->CurrentImageMapper->GetInput());
    }
  if (data)
    {
    // Compute the distance from the origin to the slice
    double origin[3];
    data->GetOrigin(origin);
    double d = vtkMath::Dot(origin, normal) + normal[3];
    // Compute the spacing to use
    double spacing[3];
    data->GetSpacing(spacing);
    double wx = normal[0]*normal[0];
    double wy = normal[1]*normal[1];
    double wz = normal[2]*normal[2];
    double s = fabs(spacing[0])*wx + fabs(spacing[1])*wy + fabs(spacing[2])*wz;
    // Round to get the nearest slice index
    int n = vtkMath::Floor(d/s + 0.5);
    n += delta;
    // Apply some limits
    int extent[6];
    data->GetWholeExtent(extent);
    int lo = vtkMath::Floor(extent[0]*wx + extent[2]*wy + extent[4]*wz + 0.5);
    int hi = vtkMath::Floor(extent[1]*wx + extent[3]*wy + extent[5]*wz + 0.5);
    n = (n < lo ? lo : n);
    n = (n > hi ? hi : n);
    // Adjust the plane
    normal[3] -= n*s - d;
    }

  // Convert the plane back to world coordinates
  if (this->CurrentImageMatrix)
    {
    double matrix[16];
    vtkMatrix4x4::Invert(*this->CurrentImageMatrix->Element, matrix);
    vtkMatrix4x4::Transpose(matrix, matrix);
    vtkMatrix4x4::MultiplyPoint(matrix, normal, normal);
    }

  // project the focal point onto this new plane
  double f = vtkMath::Dot(normal, point) + normal[3];
  point[0] += f*normal[0];
  point[1] += f*normal[1];
  point[2] += f*normal[2];

  camera->SetFocalPoint(point);
}
