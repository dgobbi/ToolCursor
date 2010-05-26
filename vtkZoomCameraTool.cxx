/*=========================================================================

  Program:   ToolCursor
  Module:    vtkZoomCameraTool.cxx

  Copyright (c) 2010 David Gobbi
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkZoomCameraTool.h"
#include "vtkObjectFactory.h"

#include "vtkToolCursor.h"
#include "vtkCamera.h"
#include "vtkRenderer.h"
#include "vtkTransform.h"
#include "vtkMath.h"

#include "vtkVolumePicker.h"

vtkStandardNewMacro(vtkZoomCameraTool);

//----------------------------------------------------------------------------
vtkZoomCameraTool::vtkZoomCameraTool()
{
  this->ZoomByDolly = 1;

  this->Transform = vtkTransform::New();
}

//----------------------------------------------------------------------------
vtkZoomCameraTool::~vtkZoomCameraTool()
{
  this->Transform->Delete();
}

//----------------------------------------------------------------------------
void vtkZoomCameraTool::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "ZoomByDolly: " << (this->ZoomByDolly ? "On\n" : "Off\n");
}

//----------------------------------------------------------------------------
void vtkZoomCameraTool::StartAction()
{
  this->Superclass::StartAction();

  vtkToolCursor *cursor = this->GetToolCursor();
  vtkCamera *camera = cursor->GetRenderer()->GetActiveCamera();

  camera->GetPosition(this->StartCameraPosition);
  camera->GetClippingRange(this->StartClippingRange);
  this->StartParallelScale = camera->GetParallelScale();
  this->StartViewAngle = camera->GetViewAngle();

  this->ZoomFactor = 1.0;

  this->Transform->Identity();
}

//----------------------------------------------------------------------------
void vtkZoomCameraTool::StopAction()
{
  this->Superclass::StopAction();
}

//----------------------------------------------------------------------------
void vtkZoomCameraTool::DoAction()
{
  this->Superclass::DoAction();

  vtkToolCursor *cursor = this->GetToolCursor();
  vtkCamera *camera = cursor->GetRenderer()->GetActiveCamera();
  vtkMatrix4x4 *viewMatrix = camera->GetViewTransformMatrix();

  // Get the camera's z axis
  double cvz[3];
  for (int i = 0; i < 3; i++)
    {
    cvz[i] = viewMatrix->GetElement(2, i);
    }

  double f[3];
  camera->GetFocalPoint(f);

  double c[3];
  c[0] = this->StartCameraPosition[0];
  c[1] = this->StartCameraPosition[1];
  c[2] = this->StartCameraPosition[2];

  // Get the initial point.
  double p0[3];
  this->GetStartPosition(p0);

  // Get the depth.
  double x, y, z;
  this->WorldToDisplay(p0, x, y, z);

  // Get the display position.
  double p[3];
  this->GetDisplayPosition(x, y);
  this->DisplayToWorld(x, y, z, p);

  // Find positions relative to camera position.
  double u[3];
  u[0] = p0[0] - f[0];
  u[1] = p0[1] - f[1];
  u[2] = p0[2] - f[2];

  // Distance from focal plane
  double df = vtkMath::Dot(u, cvz);

  // Point about which magnification occurs
  double g[3];
  g[0] = f[0] + df*cvz[0];
  g[1] = f[1] + df*cvz[1];
  g[2] = f[2] + df*cvz[2];

  // Distance from center for the two points
  double r1 = sqrt(vtkMath::Distance2BetweenPoints(g, p0));
  double r2 = sqrt(vtkMath::Distance2BetweenPoints(g, p));

  // Get the camera position
  camera->GetPosition(p);
  double dp = sqrt(vtkMath::Distance2BetweenPoints(p,g));

  // Get viewport height at the current depth
  double height = 1;
  if (camera->GetParallelProjection())
    {
    height = camera->GetParallelScale();
    }
  else
    {
    double angle = vtkMath::RadiansFromDegrees(camera->GetViewAngle());
    height = 2*dp*sin(angle/2);
    }

  // Constrain the values when they are close to the center, in order to
  // avoid magifications of zero or infinity
  double halfpi = 0.5*vtkMath::DoublePi();
  double r0 = 0.1*height;
  if (r1 < r0)
    {
    r1 = r0*(1.0 - sin((1.0 - r1/r0)*halfpi)/halfpi);
    }
  if (r2 < r0)
    {
    r2 = r0*(1.0 - sin((1.0 - r2/r0)*halfpi)/halfpi);
    }

  // Compute magnification and corresponding camera motion
  double mag = r2/r1;
  double delta = dp - dp/mag;

  this->Transform->PostMultiply();
  this->Transform->Translate(-delta*cvz[0], -delta*cvz[1], -delta*cvz[2]);

  this->ZoomFactor *= mag;

  if (camera->GetParallelProjection())
    {
    camera->SetParallelScale(this->StartParallelScale/this->ZoomFactor);
    }
  else
    {
    if (this->ZoomByDolly)
      {
      double cameraPos[3];
      this->Transform->TransformPoint(this->StartCameraPosition, cameraPos);

      camera->SetPosition(cameraPos);

      double v[3];
      v[0] = cameraPos[0] - this->StartCameraPosition[0];
      v[1] = cameraPos[1] - this->StartCameraPosition[1];
      v[2] = cameraPos[2] - this->StartCameraPosition[2];

      double dist = vtkMath::Dot(v, cvz);

      double d1 = this->StartClippingRange[0] + dist;
      double d2 = this->StartClippingRange[1] + dist;

      double tol = cursor->GetRenderer()->GetNearClippingPlaneTolerance();

      if (d1 < d2*tol)
        {
        d1 = d2*tol;
        }

      camera->SetClippingRange(d1, d2);
      }
    else
      {
      // Zoom by changing the view angle

      double h = 2*tan(0.5*vtkMath::RadiansFromDegrees(this->StartViewAngle));
      h /= this->ZoomFactor;
      double viewAngle = 2*vtkMath::DegreesFromRadians(atan(0.5*h));

      camera->SetViewAngle(viewAngle);
      }
    }
}

//----------------------------------------------------------------------------
void vtkZoomCameraTool::ConstrainCursor(double position[3],
                                          double normal[3])
{
  vtkToolCursor *cursor = this->GetToolCursor();
  vtkCamera *camera = cursor->GetRenderer()->GetActiveCamera();
  vtkMatrix4x4 *matrix = camera->GetViewTransformMatrix();

  // Force the normal to point towards camera
  normal[0] = matrix->GetElement(2, 0);
  normal[1] = matrix->GetElement(2, 1);
  normal[2] = matrix->GetElement(2, 2);
}

