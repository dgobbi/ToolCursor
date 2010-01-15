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

vtkCxxRevisionMacro(vtkRotateCameraAction, "$Revision: 1.4 $");
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

  // Get the camera
  vtkSurfaceCursor *cursor = this->GetSurfaceCursor();
  vtkCamera *camera = cursor->GetRenderer()->GetActiveCamera();
  vtkMatrix4x4 *viewMatrix = camera->GetViewTransformMatrix();

  // Get the camera's x, y, and z axes
  double cvx[3], cvy[3], cvz[3];
  for (int i = 0; i < 3; i++)
    {
    cvx[i] = viewMatrix->GetElement(0, i);
    cvy[i] = viewMatrix->GetElement(1, i);
    cvz[i] = viewMatrix->GetElement(2, i);
    }

  // Get the camera's position
  double cameraPos[3];
  camera->GetPosition(cameraPos);

  // Center of rotation (focal point) and rotation sphere radius
  double f[3];
  this->GetCenterOfRotation(f);
  double r = this->GetRadius();

  // Camera distance, squared, and rotation radius, squared
  double d2 = vtkMath::Distance2BetweenPoints(cameraPos,f);
  double r2 = r*r;

  // The interaction uses a plane parallel to the focal plane
  // but far enough in front of the focal plane that it defines
  // a circle that contains 75% of the visible sphere area.
  double s2max = 0.75*r2/d2*(d2 - r2);

  // And here is the position of the desired interaction plane.
  double fg = sqrt(r2 - s2max);
  double g[3];
  g[0] = f[0] + fg*cvz[0];
  g[1] = f[1] + fg*cvz[1];
  g[2] = f[2] + fg*cvz[2];

  // Get the current display position. 
  double x, y;
  this->GetDisplayPosition(x, y);

  // Get the view ray and see where it intersects the sphere of rotation.
  // This involves several steps and solving a quadratic equation.
  double p1[3], p2[3];
  this->DisplayToWorld(x, y, 0.0, p1);
  this->DisplayToWorld(x, y, 1.0, p2);

  // Vector along direction of view-ray line
  double v[3];
  v[0] = p2[0] - p1[0];
  v[1] = p2[1] - p1[1];
  v[2] = p2[2] - p1[2];

  // Vector from center of rotation to first line point
  double u[3];
  u[0] = p1[0] - f[0];
  u[1] = p1[1] - f[1];
  u[2] = p1[2] - f[2];

  // Here are the coefficients of the quadratic equation, at^2 + bt + c = 0,
  // which will give us the parametric distance "t" along the view ray line
  double a = vtkMath::Dot(v, v);
  double b = 2*vtkMath::Dot(u, v);
  double c = vtkMath::Dot(u, u) - r*r;

  // The value under the square root in the solution to the quadratic
  double discriminant = b*b - 4*a*c;

  double t = 1.0;
  if (discriminant > 0)
    {
    // The discriminant is positive, so there are two real solutions.
    // Take the smaller of the two roots.
    t = (-b - sqrt(discriminant))/(2*a);
    }
  //else { cerr << "outside 1\n"; }

  // Drop down to intersect with the interaction plane.
  u[0] = p1[0] - g[0];
  u[1] = p1[1] - g[1];
  u[2] = p1[2] - g[2];
  double tplane = -vtkMath::Dot(u,cvz)/vtkMath::Dot(v,cvz);
  if (t > tplane)
    {
    //cerr << "in the ring 1\n";
    t = tplane;
    }

  // Compute the intersection point
  double p[3];
  p[0] = p1[0]*(1 - t) + p2[0]*t;
  p[1] = p1[1]*(1 - t) + p2[1]*t;
  p[2] = p1[2]*(1 - t) + p2[2]*t;

  // Displacement from center of interaction plane
  double vi[3];
  vi[0] = p[0] - g[0];
  vi[1] = p[1] - g[1];
  vi[2] = p[2] - g[2];

  // Camera coords of point, centered at center of rotation
  double cx = vtkMath::Dot(vi, cvx);
  double cy = vtkMath::Dot(vi, cvy);
  double s2 = cx*cx + cy*cy;

  double w2[3];
  w2[0] = cvz[0];
  w2[1] = cvz[1];
  w2[2] = cvz[2];

  // If point is outside the circle that defines 75% of the visible sphere
  //cerr << "s = " << sqrt(s2) << ", smax = " << sqrt(s2max) << "\n"; 
  if (s2 > s2max)
    {
    // Get the current display position. 
    double x, y;
    this->GetLastDisplayPosition(x, y);

    // Get the view ray and see where it intersects the sphere of rotation.
    // This involves several steps and solving a quadratic equation.
    double p1[3], p2[3];
    this->DisplayToWorld(x, y, 0.0, p1);
    this->DisplayToWorld(x, y, 1.0, p2);

    // Compute the intersection point
    double pl[3];
    pl[0] = p1[0]*(1 - t) + p2[0]*t;
    pl[1] = p1[1]*(1 - t) + p2[1]*t;
    pl[2] = p1[2]*(1 - t) + p2[2]*t;

    //if (s2 > r2/d2*(d2 - r2)) { cerr << "outside 2\n"; }
    //else { cerr << "in the ring 2\n"; }
    double s = sqrt(s2);
    double smax = sqrt(s2max);

    double cosphi = cx/s;
    double sinphi = cy/s;

    double psi = 2*(s - smax)/r + asin(smax/r);

    //cerr << "psi = " << psi << ", phi = " << atan2(sinphi,cosphi) << "\n";
    double q = r*sin(psi);
    double sx = q*cosphi;
    double sy = q*sinphi;
    double sz = r*cos(psi);

    p[0] = sx*cvx[0] + sy*cvy[0] + sz*cvz[0] + f[0];
    p[1] = sx*cvx[1] + sy*cvy[1] + sz*cvz[1] + f[1];
    p[2] = sx*cvx[2] + sy*cvy[2] + sz*cvz[2] + f[2];
    }

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
  double delta = vtkMath::Norm(w);
  if (delta/r < 1e-7)
    {
    //cerr << "same point!\n";
    return;
    }

  // The line between pr and our center-of-rotation is the line of rotation
  double pr[3];

  // First, get the direction of the line of rotation
  double n[3];
  vtkMath::Cross(w, w2, n);
  double nlen = vtkMath::Norm(n);

  if (nlen/delta < 1e-5 || s2 > s2max)
    {
    // Can't define a line
    pr[0] = f[0];
    pr[1] = f[1];
    pr[2] = f[2];
    }
  else
    {
    n[0] = n[0] / nlen;
    n[1] = n[1] / nlen;
    n[2] = n[2] / nlen;

    // Find the point on the line that is closest to the two points.
    u[0] = p0[0] - f[0];
    u[1] = p0[1] - f[1];
    u[2] = p0[2] - f[2];

    t = vtkMath::Dot(u, n);
    pr[0] = f[0] + t*n[0];
    pr[1] = f[1] + t*n[1];
    pr[2] = f[2] + t*n[2];
    }

  // The point pr and the points p, p0 form our rotation angle
  double v1[3];
  v1[0] = p0[0] - pr[0];
  v1[1] = p0[1] - pr[1];
  v1[2] = p0[2] - pr[2];

  double v2[3];
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

  this->Transform->TransformPoint(this->StartCameraPosition, cameraPos);
  this->Transform->TransformVector(this->StartCameraViewUp, cvy);

  camera->SetPosition(cameraPos);
  camera->SetFocalPoint(f);
  camera->SetViewUp(cvy);
}

//----------------------------------------------------------------------------
void vtkRotateCameraAction::ConstrainCursor(double *position,
                                            double normal[3]) 
{
  // Use the original normal.
  this->GetNormal(normal);
}

