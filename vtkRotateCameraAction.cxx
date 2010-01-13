/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkRotateCameraAction.cxx,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkRotateCameraAction.h"
#include "vtkObjectFactory.h"

#include "vtkSurfaceCursor.h"
#include "vtkCamera.h"
#include "vtkRenderer.h"
#include "vtkTransform.h"
#include "vtkMath.h"

#include "vtkVolumePicker.h"

vtkCxxRevisionMacro(vtkRotateCameraAction, "$Revision: 1.1 $");
vtkStandardNewMacro(vtkRotateCameraAction);

//----------------------------------------------------------------------------
vtkRotateCameraAction::vtkRotateCameraAction()
{
  this->Transform = vtkTransform::New();
}

//----------------------------------------------------------------------------
vtkRotateCameraAction::~vtkRotateCameraAction()
{
  this->Transform->Delete();
}

//----------------------------------------------------------------------------
void vtkRotateCameraAction::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

//----------------------------------------------------------------------------
void vtkRotateCameraAction::StartAction()
{
  this->Superclass::StartAction();

  vtkSurfaceCursor *cursor = this->GetSurfaceCursor();
  vtkCamera *camera = cursor->GetRenderer()->GetActiveCamera();

  cursor->GetNormal(this->Normal);
  camera->GetFocalPoint(this->CenterOfRotation);
  camera->GetPosition(this->StartCameraPosition);
  camera->GetViewUp(this->StartCameraViewUp);

  this->Radius = sqrt(vtkMath::Distance2BetweenPoints(
    this->CenterOfRotation, this->GetStartPosition()));

  this->Transform->Identity();
} 

//----------------------------------------------------------------------------
void vtkRotateCameraAction::StopAction()
{
  this->Superclass::StopAction();
}

//----------------------------------------------------------------------------
void vtkRotateCameraAction::DoAction()
{
  this->Superclass::DoAction();

  // Get the current display position. 
  double x, y;
  this->GetDisplayPosition(x, y);

  // Get the view ray and see where it intersects the sphere of rotation.
  double p1[3], p2[3];
  this->DisplayToWorld(x, y, 0.0, p1);
  this->DisplayToWorld(x, y, 1.0, p2);

  // Center of rotation and radius
  double f[3];
  this->GetCenterOfRotation(f);
  double r = this->GetRadius();

  // Vector from center of rotation to first line point
  double u[3];
  u[0] = p1[0] - f[0];
  u[1] = p1[1] - f[1];
  u[2] = p1[2] - f[2];

  // Vector along direction of line
  double v[3];
  v[0] = p2[0] - p1[0];
  v[1] = p2[1] - p1[1];
  v[2] = p2[2] - p1[2];

  // Here are the coefficients of the quadratic equation, at^2 + bt + c = 0,
  // which will give us the parametric distance "t" along the view ray line
  double a = vtkMath::Dot(v, v);
  double b = 2*vtkMath::Dot(u, v);
  double c = vtkMath::Dot(u, u) - r*r;

  // The value under the square root in the solution to the quadratic
  double d = b*b - 4*a*c;

  double t = 0;
  if (d > 0)
    {
    // Take the smaller of the two roots.
    t = (-b - sqrt(d))/(2*a);
    }
  else
    {
    // The position is off the sphere
    t = -b/(2*a);
    }

  // Finally, compute the intersection point
  double p[3];
  p[0] = p1[0]*(1 - t) + p2[0]*t;
  p[1] = p1[1]*(1 - t) + p2[1]*t;
  p[2] = p1[2]*(1 - t) + p2[2]*t;

  // Get the initial point.
  double p0[3];
  this->GetStartPosition(p0);

  //cerr << " p0 (" << p0[0] << " " << p0[1] << " " << p0[2] << ")";
  //cerr << " t " << t << " ";
  //cerr << " p (" << p[0] << " " << p[1] << " " << p[2] << ")\n";

  // Okay, now we have our "Start" point, and our new point.
  // Find the vector between them.
  double w[3];
  w[0] = p[0] - p0[0];
  w[1] = p[1] - p0[1];
  w[2] = p[2] - p0[2];

  if (vtkMath::Dot(w, w) < 1e-20)
    {
    //cerr << "same point!\n";
    return;
    }

  // Cross with the view ray to get the desired rotation axis.
  double n[3];
  vtkMath::Cross(w, v, n);
  vtkMath::Normalize(n);

  // The rotation axis, together with the camera focal point, form a line.
  // Find the point on that line that is closest to the two points.
  u[0] = p0[0] - f[0]; 
  u[1] = p0[1] - f[1]; 
  u[2] = p0[2] - f[2]; 

  double pr[3];
  double s = vtkMath::Dot(u, n);
  pr[0] = f[0] + s*n[0];
  pr[1] = f[1] + s*n[1];
  pr[2] = f[2] + s*n[2];

  // The point pr and the points p, p0 form our rotation angle
  double v1[3], v2[3];
  v1[0] = p0[0] - pr[0];
  v1[1] = p0[1] - pr[1];
  v1[2] = p0[2] - pr[2];
  v2[0] = p[0] - pr[0];
  v2[1] = p[1] - pr[1];
  v2[2] = p[2] - pr[2];

  double v1n = vtkMath::Norm(v1);
  double v2n = vtkMath::Norm(v2);

  double v3[3];
  vtkMath::Cross(v1, v2, v3);
  double v3n = vtkMath::Norm(v3);

  //cerr << "v1n = " << v1n << ", v2n = " << v2n << ", v3n = " << v3n << "\n";

  double sintheta = v3n/(v1n*v2n);
  double costheta = vtkMath::Dot(v1, v2)/(v1n*v2n);

  // Finally: the angle and the rotation axis
  double angle = atan2(sintheta, costheta);
  v3[0] /= v3n;
  v3[1] /= v3n;
  v3[2] /= v3n;

  //cerr << "sin " << sintheta << ", cos " << costheta << ", angle " << vtkMath::DegreesFromRadians(angle) << "\n";

  //cerr << "n v3 (" << n[0] << " " << v3[0] << ", " << n[1] << " " << v3[1] << ", " << n[2] << " " << v3[2] << ")\n";

  // Rotate the original direction of projection and view-up
  this->Transform->PostMultiply();
  this->Transform->Translate(-f[0], -f[1], -f[2]);
  this->Transform->RotateWXYZ(vtkMath::DegreesFromRadians(-angle), v3);
  this->Transform->Translate(f[0], f[1], f[2]);

  vtkSurfaceCursor *cursor = this->GetSurfaceCursor();
  vtkCamera *camera = cursor->GetRenderer()->GetActiveCamera();

  double cameraPos[3], cameraViewUp[3];
  this->Transform->TransformPoint(this->StartCameraPosition, cameraPos);
  this->Transform->TransformVector(this->StartCameraViewUp, cameraViewUp);

  camera->SetPosition(cameraPos);
  camera->SetFocalPoint(f);
  camera->SetViewUp(cameraViewUp);
}

//----------------------------------------------------------------------------
void vtkRotateCameraAction::ConstrainCursor(double *position,
                                            double normal[3]) 
{
  // Use the original normal.
  this->GetNormal(normal);
}

