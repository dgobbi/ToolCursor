/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkPushPlaneAction.cxx,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkPushPlaneAction.h"
#include "vtkObjectFactory.h"

#include "vtkSurfaceCursor.h"
#include "vtkRenderer.h"
#include "vtkVolumePicker.h"
#include "vtkProp3DCollection.h"
#include "vtkImageActor.h"
#include "vtkVolumeMapper.h"

vtkCxxRevisionMacro(vtkPushPlaneAction, "$Revision: 1.1 $");
vtkStandardNewMacro(vtkPushPlaneAction);

//----------------------------------------------------------------------------
vtkPushPlaneAction::vtkPushPlaneAction()
{
  this->ImageActor = 0;
  this->VolumeMapper = 0;
  this->Mapper = 0;
  this->PlaneId = -1;
}

//----------------------------------------------------------------------------
vtkPushPlaneAction::~vtkPushPlaneAction()
{
}

//----------------------------------------------------------------------------
void vtkPushPlaneAction::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

//----------------------------------------------------------------------------
void vtkPushPlaneAction::StartAction()
{
  this->Superclass::StartAction();

  // Get all the object needed for the interaction
  vtkVolumePicker *picker = this->SurfaceCursor->GetPicker();
  vtkProp3DCollection *props = picker->GetProp3Ds();
  vtkCollectionSimpleIterator pit;
  props->InitTraversal(pit);
  this->ImageActor = vtkImageActor::SafeDownCast(props->GetNextProp3D(pit));
  this->Mapper = picker->GetMapper();
  this->VolumeMapper = vtkVolumeMapper::SafeDownCast(this->Mapper);
  if (this->VolumeMapper && picker->GetCroppingPlaneId() >= 0)
    {
    this->Mapper = 0;
    this->PlaneId = picker->GetCroppingPlaneId();
    }
  else
    {
    this->VolumeMapper = 0;
    this->PlaneId = picker->GetClippingPlaneId();
    }
} 

//----------------------------------------------------------------------------
void vtkPushPlaneAction::StopAction()
{
  this->Superclass::StopAction();
}

//----------------------------------------------------------------------------
void vtkPushPlaneAction::DoAction()
{
  this->Superclass::DoAction();

  if (this->ImageActor)
    {
    }
  else if (this->Mapper)
    {
    }
  else if (this->VolumeMapper)
    {
    cerr << "Cropping Plane " << this->PlaneId << "\n";
    }
}

