/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkPanCameraAction.cxx,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkPanCameraAction.h"
#include "vtkObjectFactory.h"

#include "vtkSurfaceCursor.h"
#include "vtkCamera.h"
#include "vtkRenderer.h"
#include "vtkTransform.h"
#include "vtkMath.h"

#include "vtkVolumePicker.h"

vtkCxxRevisionMacro(vtkPanCameraAction, "$Revision: 1.1 $");
vtkStandardNewMacro(vtkPanCameraAction);

//----------------------------------------------------------------------------
vtkPanCameraAction::vtkPanCameraAction()
{
  this->Transform = vtkTransform::New();
}

//----------------------------------------------------------------------------
vtkPanCameraAction::~vtkPanCameraAction()
{
  this->Transform->Delete();
}

//----------------------------------------------------------------------------
void vtkPanCameraAction::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

//----------------------------------------------------------------------------
void vtkPanCameraAction::StartAction()
{
  this->Superclass::StartAction();

  vtkSurfaceCursor *cursor = this->GetSurfaceCursor();
  vtkCamera *camera = cursor->GetRenderer()->GetActiveCamera();

  camera->GetFocalPoint(this->StartCameraFocalPoint);
  camera->GetPosition(this->StartCameraPosition);

  this->Transform->Identity();
} 

//----------------------------------------------------------------------------
void vtkPanCameraAction::StopAction()
{
  this->Superclass::StopAction();
}

//----------------------------------------------------------------------------
void vtkPanCameraAction::DoAction()
{
  this->Superclass::DoAction();

  vtkSurfaceCursor *cursor = this->GetSurfaceCursor();
  vtkCamera *camera = cursor->GetRenderer()->GetActiveCamera();

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

  // Get the vector.
  double v[3];
  v[0] = p[0] - p0[0];
  v[1] = p[1] - p0[1];
  v[2] = p[2] - p0[2];

  this->Transform->PostMultiply();
  this->Transform->Translate(-v[0], -v[1], -v[2]);

  double cameraPos[3], cameraFocalPoint[3];
  this->Transform->TransformPoint(this->StartCameraPosition, cameraPos);
  this->Transform->TransformPoint(this->StartCameraFocalPoint,
                                  cameraFocalPoint);

  camera->SetPosition(cameraPos);
  camera->SetFocalPoint(cameraFocalPoint);
}

