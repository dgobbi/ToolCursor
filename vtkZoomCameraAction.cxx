/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkZoomCameraAction.cxx,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkZoomCameraAction.h"
#include "vtkObjectFactory.h"

#include "vtkSurfaceCursor.h"
#include "vtkCamera.h"
#include "vtkRenderer.h"
#include "vtkTransform.h"
#include "vtkMath.h"

#include "vtkVolumePicker.h"

vtkCxxRevisionMacro(vtkZoomCameraAction, "$Revision: 1.1 $");
vtkStandardNewMacro(vtkZoomCameraAction);

//----------------------------------------------------------------------------
vtkZoomCameraAction::vtkZoomCameraAction()
{
  this->Transform = vtkTransform::New();
}

//----------------------------------------------------------------------------
vtkZoomCameraAction::~vtkZoomCameraAction()
{
  this->Transform->Delete();
}

//----------------------------------------------------------------------------
void vtkZoomCameraAction::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

//----------------------------------------------------------------------------
void vtkZoomCameraAction::StartAction()
{
  this->Superclass::StartAction();

  vtkSurfaceCursor *cursor = this->GetSurfaceCursor();
  vtkCamera *camera = cursor->GetRenderer()->GetActiveCamera();

  camera->GetPosition(this->StartCameraPosition);
  camera->GetClippingRange(this->StartClippingRange);
  this->StartParallelScale = camera->GetParallelScale();
  this->ZoomFactor = 1.0;

  this->Transform->Identity();
} 

//----------------------------------------------------------------------------
void vtkZoomCameraAction::StopAction()
{
  this->Superclass::StopAction();
}

//----------------------------------------------------------------------------
void vtkZoomCameraAction::DoAction()
{
  this->Superclass::DoAction();

  vtkSurfaceCursor *cursor = this->GetSurfaceCursor();
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

  // Magnification factor
  double mag = sqrt(vtkMath::Distance2BetweenPoints(g, p)/
                    vtkMath::Distance2BetweenPoints(g, p0));

  // Get the camera position to determine necessary motion.
  camera->GetPosition(p);
  double dp = sqrt(vtkMath::Distance2BetweenPoints(p,g));
  double delta = dp - dp/mag;

  this->Transform->PostMultiply();
  this->Transform->Translate(-delta*cvz[0], -delta*cvz[1], -delta*cvz[2]);

  this->ZoomFactor *= mag;

  if (camera->GetParallelProjection())
    {
    camera->SetParallelScale(this->StartParallelScale*this->ZoomFactor);
    }
  else
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
}

