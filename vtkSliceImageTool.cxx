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
    double d = 0.0;
    if ((cursor->GetModifier() & VTK_TOOL_WHEEL_BWD) != 0)
      {
      d = -1.0;
      }
    else if ((cursor->GetModifier() & VTK_TOOL_WHEEL_FWD) != 0)
      {
      d = +1.0;
      }
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

