/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkSurfaceCursorAction.cxx,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkSurfaceCursorAction.h"
#include "vtkObjectFactory.h"

#include "vtkSurfaceCursor.h"
#include "vtkRenderer.h"
#include "vtkCamera.h"
#include "vtkVolumePicker.h"
#include "vtkMath.h"

vtkCxxRevisionMacro(vtkSurfaceCursorAction, "$Revision: 1.2 $");
vtkStandardNewMacro(vtkSurfaceCursorAction);

//----------------------------------------------------------------------------
vtkSurfaceCursorAction::vtkSurfaceCursorAction()
{
  this->SurfaceCursor = 0;
}

//----------------------------------------------------------------------------
vtkSurfaceCursorAction::~vtkSurfaceCursorAction()
{
  // SurfaceCursor is not reference counted and therefore not deleted.
  this->SurfaceCursor = 0;
}

//----------------------------------------------------------------------------
void vtkSurfaceCursorAction::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

//----------------------------------------------------------------------------
void vtkSurfaceCursorAction::SetSurfaceCursor(vtkSurfaceCursor *cursor)
{
  if (cursor != this->SurfaceCursor)
    {
    this->SurfaceCursor = cursor;
    this->Modified();
    }
}

//----------------------------------------------------------------------------
void vtkSurfaceCursorAction::StartAction()
{
  this->SurfaceCursor->GetPosition(this->StartPosition);
  this->SurfaceCursor->GetDisplayPosition(this->StartDisplayPosition);
  this->DisplayPosition[0] = this->StartDisplayPosition[0];
  this->DisplayPosition[1] = this->StartDisplayPosition[1];
  this->LastDisplayPosition[0] = this->DisplayPosition[0];
  this->LastDisplayPosition[1] = this->DisplayPosition[1];
} 

//----------------------------------------------------------------------------
void vtkSurfaceCursorAction::StopAction()
{
}

//----------------------------------------------------------------------------
void vtkSurfaceCursorAction::DoAction()
{
  this->LastDisplayPosition[0] = this->DisplayPosition[0];
  this->LastDisplayPosition[1] = this->DisplayPosition[1];
  this->SurfaceCursor->GetDisplayPosition(this->DisplayPosition);
}

//----------------------------------------------------------------------------
void vtkSurfaceCursorAction::WorldToDisplay(const double world[3],
                                            double &x, double &y, double &z)
{
  vtkRenderer *renderer = this->SurfaceCursor->GetRenderer();
  double hcoord[4];
  hcoord[0] = world[0];
  hcoord[1] = world[1];
  hcoord[2] = world[2];
  hcoord[3] = 1.0;

  // Use the horrendous viewport interterface for conversions.
  renderer->SetWorldPoint(hcoord);
  renderer->WorldToDisplay();
  renderer->GetDisplayPoint(hcoord);
  x = hcoord[0];
  y = hcoord[1];
  z = hcoord[2];
}

//----------------------------------------------------------------------------
void vtkSurfaceCursorAction::DisplayToWorld(double x, double y, double z,
                                            double world[3])
{
  vtkRenderer *renderer = this->SurfaceCursor->GetRenderer();
  double hcoord[4];

  // Use the viewport interterface for conversions.
  renderer->SetDisplayPoint(x, y, z);
  renderer->DisplayToWorld();
  renderer->GetWorldPoint(hcoord);
  world[0] = hcoord[0]/hcoord[3];
  world[1] = hcoord[1]/hcoord[3];
  world[2] = hcoord[2]/hcoord[3];
}

//----------------------------------------------------------------------------
void vtkSurfaceCursorAction::GetViewRay(double x, double y, double z,
                                        double p[3], double v[3])
{
  double p1[3], p2[3];
  this->DisplayToWorld(x, y, 0.0, p1);
  this->DisplayToWorld(x, y, z, p);
  this->DisplayToWorld(x, y, 1.0, p2);

  v[0] = p2[0] - p1[0];
  v[1] = p2[1] - p1[1];
  v[2] = p2[2] - p1[2];

  vtkMath::Normalize(v);
}

