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
#include "vtkVolumePicker.h"


vtkCxxRevisionMacro(vtkSurfaceCursorAction, "$Revision: 1.1 $");
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

