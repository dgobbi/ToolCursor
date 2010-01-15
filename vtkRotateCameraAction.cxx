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

vtkCxxRevisionMacro(vtkRotateCameraAction, "$Revision: 1.5 $");
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

  camera->GetFocalPoint(this->CenterOfRotation);
  camera->GetPosition(this->StartCameraPosition);
  camera->GetViewUp(this->StartCameraViewUp);

  // The radius is the initial distance from the focal point
  this->Radius = sqrt(vtkMath::Distance2BetweenPoints(
    this->CenterOfRotation, this->GetStartPosition()));

  // Pre-compute the view plane normal and the camera distance
  double *f = this->CenterOfRotation;
  double *p = this->StartCameraPosition;

  double cvz[3];
  cvz[0] = p[0] - f[0];
  cvz[1] = p[1] - f[1];
  cvz[2] = p[2] - f[2];

  double d = vtkMath::Norm(cvz);
  cvz[0] /= d;
  cvz[1] /= d;
  cvz[2] /= d;

  // "Sticky" means stick the cursor to a sphere of radius r
  this->Sticky = 1;

  // We stick the cursor to the start position
  this->GetStartPosition(this->StickyPosition);

  // Check for conditions where "sticky" interaction won't work

  // Sticky depends on there being an object to stick to
  if (cursor->GetPickFlags() == 0)
    {
    this->Sticky = 0;
    }

  // Get viewport height at focal plane
  double height = 1;
  if (camera->GetParallelProjection())
    {
    height = camera->GetParallelScale();
    }
  else
    {
    height = 2*d*sin(vtkMath::RadiansFromDegrees(camera->GetViewAngle()/2));
    }

  if (this->Radius < height*0.05)
    {
    this->Radius = height*0.05;
    this->Sticky = 0;
    }
  else
    {
    // Can only do "sticky" if initial point is in front of focal plane

    // Get the distance of the point from the focal plane
    double p0[3];
    this->GetStartPosition(p0);

    double v[3];
    v[0] = p0[0] - f[0];
    v[1] = p0[1] - f[1];
    v[2] = p0[2] - f[2];

    double df = vtkMath::Dot(v, cvz);

    // The "r/d" check for perspective projections is needed because
    // the edge of the "rotation sphere" from the camera's viewpoint
    // will actually be in front of the center of the sphere.
    if ((camera->GetParallelProjection() && df < this->Radius/d) ||
        df < 0)
      {
      this->Sticky = 0;
      }
    }

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

  // Here are the values that we will be setting in the code: the angle,
  // in radians, and the axis of rotation
  double rotationAngle = 0;
  double rotationAxis[3];
  rotationAxis[0] = 0;
  rotationAxis[1] = 0;
  rotationAxis[2] = 1;

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

  // Camera distance, squared, and rotation radius, squared.
  double d2 = vtkMath::Distance2BetweenPoints(cameraPos,f);
  double r2 = r*r;

  // The interaction uses a plane parallel to the focal plane
  // but far enough in front of the focal plane that it defines
  // a circle that contains a certain percentage of the visible
  // sphere area, say 99%.
  double s2max = 0.99*r2/d2*(d2 - r2);

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

  // Vector to the center of interaction plane.
  double u[3];
  u[0] = p1[0] - g[0];
  u[1] = p1[1] - g[1];
  u[2] = p1[2] - g[2];

  // Find the "t" value for the interaction plane
  double t = -vtkMath::Dot(u,cvz)/vtkMath::Dot(v,cvz);

  if (this->Sticky)
    {
    // Vector from center of rotation to first line point
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

    if (discriminant > 0)
      {
      // The discriminant is positive, so there are two real solutions.
      // Take the smaller of the two roots.
      double tSphere = (-b - sqrt(discriminant))/(2*a);

      // Only allow "t" values that are in front of interaction plane
      if (tSphere < t)
        {
        t = tSphere;
        }
      }

    // Use "t" to compute the intersection point for the view ray
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

    // If point is off the sphere and on the interaction plane instead
    if (s2 > s2max)
      {
      // Map the off-sphere point to a position on the sphere
      double s = sqrt(s2);
      double smax = sqrt(s2max);

      double cosphi = cx/s;
      double sinphi = cy/s;

      // Convert "s" to an angle of lattitude on the sphere, in a way
      // such that lattitude values are continuous when moving onto and
      // off of the sphere
      double psi = 2*(s - smax)/r + asin(smax/r);

      // Use spherical coordinates to compute position on sphere
      double q = r*sin(psi);
      double sx = q*cosphi;
      double sy = q*sinphi;
      double sz = r*cos(psi);

      // Apply camera transform to put into world coordinates
      p[0] = sx*cvx[0] + sy*cvy[0] + sz*cvz[0] + f[0];
      p[1] = sx*cvx[1] + sy*cvy[1] + sz*cvz[1] + f[1];
      p[2] = sx*cvx[2] + sy*cvy[2] + sz*cvz[2] + f[2];
      }

    // Get the initial point, i.e. the point we want to stick to.
    double p0[3];
    this->GetStickyPosition(p0);

    // Find the vector between the start point and the current point
    double w[3];
    w[0] = p[0] - p0[0];
    w[1] = p[1] - p0[1];
    w[2] = p[2] - p0[2];

    // Verify that there is significant motion between the points
    double delta = vtkMath::Norm(w);
    if (delta/r < 1e-7)
      {
      return;
      }

    // pr will be the point we rotate around
    double pr[3];

    // Get the direction of the rotation vector.  It must be perpendicular
    // to the vector between the points, and for intuitive interaction it
    // must also be parallel to the view plane i.e. perpendicular to the
    // view plane normal.
    double n[3];
    vtkMath::Cross(w, cvz, n);
    double nlen = vtkMath::Norm(n);

    // If "w" and "cvz" are parallel, or if the mouse is off the rotation
    // sphere, then don't keep the rotation vector perpendicular to the
    // view plane normal.  Instead, rotate around focal point.
    if (nlen/delta < 1e-5 || s2 > s2max)
      {
      // Can't define a line
      pr[0] = f[0];
      pr[1] = f[1];
      pr[2] = f[2];
      }
    else
      {
      // If we get here, everything's fine and the point is sticky.

      // Normalize our rotation vector.
      n[0] = n[0] / nlen;
      n[1] = n[1] / nlen;
      n[2] = n[2] / nlen;

      // Project point p0 onto the rotation vector to get pr.
      u[0] = p0[0] - f[0];
      u[1] = p0[1] - f[1];
      u[2] = p0[2] - f[2];

      t = vtkMath::Dot(u, n);
      pr[0] = f[0] + t*n[0];
      pr[1] = f[1] + t*n[1];
      pr[2] = f[2] + t*n[2];
      }

    // The point pr and the points p, p0 form our rotation angle. Start
    // by computing the vectors to p and p0.
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

    // Compute the rotation angle via a cross product.
    double v3[3];
    vtkMath::Cross(v1, v2, v3);
    double v3n = vtkMath::Norm(v3);

    double sintheta = v3n/(v1n*v2n);
    double costheta = vtkMath::Dot(v1, v2)/(v1n*v2n);

    rotationAngle = atan2(sintheta, costheta);
    rotationAxis[0] = v3[0]/v3n;
    rotationAxis[1] = v3[1]/v3n;
    rotationAxis[2] = v3[2]/v3n;
    }
  else
    {
    // Here is the code for non-sticky interaction, it is much simpler.

    // Use "t" to compute the intersection point for the view ray
    double p[3];
    p[0] = p1[0]*(1 - t) + p2[0]*t;
    p[1] = p1[1]*(1 - t) + p2[1]*t;
    p[2] = p1[2]*(1 - t) + p2[2]*t;

    // Do the same for the previous position
    this->GetLastDisplayPosition(x, y);

    this->DisplayToWorld(x, y, 0.0, p1);
    this->DisplayToWorld(x, y, 1.0, p2);

    double pl[3];
    pl[0] = p1[0]*(1 - t) + p2[0]*t;
    pl[1] = p1[1]*(1 - t) + p2[1]*t;
    pl[2] = p1[2]*(1 - t) + p2[2]*t;

    // Compute the motion vector
    double w[3];
    w[0] = pl[0] - p[0];
    w[1] = pl[1] - p[1];
    w[2] = pl[2] - p[2];

    // And the rotation angle corresponding to the vector
    rotationAngle = 2*vtkMath::Norm(w)/r;

    // Rotation axis is perpendicular to view plane normal and motion vector
    vtkMath::Cross(w, cvz, rotationAxis);
    vtkMath::Normalize(rotationAxis);
    }

  // Get ready to apply the camera transformation
  this->Transform->PostMultiply();

  // First, turn the transform into a pure rotation matrix, to avoid
  // accumulation of roundoff errors in the position
  double oldTrans[3];
  this->Transform->GetPosition(oldTrans);
  this->Transform->Translate(-oldTrans[0], -oldTrans[1], -oldTrans[2]);

  // Increment by the new rotation
  this->Transform->RotateWXYZ(vtkMath::DegreesFromRadians(-rotationAngle),
                              rotationAxis);

  // Center the rotation at the focal point
  this->Transform->PreMultiply();
  this->Transform->Translate(-f[0], -f[1], -f[2]);
  this->Transform->PostMultiply();
  this->Transform->Translate(f[0], f[1], f[2]);

  // Rotate the original direction of projection and view-up
  this->Transform->TransformPoint(this->StartCameraPosition, cameraPos);
  this->Transform->TransformVector(this->StartCameraViewUp, cvy);

  camera->SetPosition(cameraPos);
  camera->SetViewUp(cvy);
}

