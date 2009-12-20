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
#include "vtkCamera.h"
#include "vtkVolumePicker.h"
#include "vtkProp3DCollection.h"
#include "vtkImageActor.h"
#include "vtkVolumeMapper.h"
#include "vtkPlaneCollection.h"
#include "vtkPlane.h"
#include "vtkTransform.h"
#include "vtkImageData.h"

vtkCxxRevisionMacro(vtkPushPlaneAction, "$Revision: 1.2 $");
vtkStandardNewMacro(vtkPushPlaneAction);

//----------------------------------------------------------------------------
vtkPushPlaneAction::vtkPushPlaneAction()
{
  this->Transform = vtkTransform::New();

  this->ImageActor = 0;
  this->VolumeMapper = 0;
  this->Mapper = 0;
  this->PlaneId = -1;
}

//----------------------------------------------------------------------------
vtkPushPlaneAction::~vtkPushPlaneAction()
{
  this->Transform->Delete();
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

  // Get all the necessary information about the picked prop.
  this->GetPropInformation();
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

  if (!this->IsPlaneValid())
    {
    return;
    }

  vtkRenderer *renderer = this->SurfaceCursor->GetRenderer();
  vtkCamera *camera = renderer->GetActiveCamera();

  double origin[3];
  double normal[3];
  this->GetOrigin(origin);
  this->GetNormal(normal);

  double cx = this->DisplayPosition[0];
  double cy = this->DisplayPosition[1];

  double px = this->LastDisplayPosition[0];
  double py = this->LastDisplayPosition[1];

  double t = px - cx;

  origin[0] = origin[0] + t*normal[0];
  origin[1] = origin[1] + t*normal[1];
  origin[2] = origin[2] + t*normal[2];

  this->SetOrigin(origin);
}

//----------------------------------------------------------------------------
void vtkPushPlaneAction::GetPropInformation()
{
  // Get all the object needed for the interaction
  vtkVolumePicker *picker = this->SurfaceCursor->GetPicker();

  vtkProp3DCollection *props = picker->GetProp3Ds();
  vtkCollectionSimpleIterator pit;
  props->InitTraversal(pit);
  vtkProp3D *prop = props->GetNextProp3D(pit);

  this->Transform->SetMatrix(prop->GetMatrix());
  this->ImageActor = vtkImageActor::SafeDownCast(prop);
  this->Mapper = picker->GetMapper();
  this->VolumeMapper = vtkVolumeMapper::SafeDownCast(this->Mapper);

  // Initialize plane to "no plane" value
  this->PlaneId = -1;

  // Is this a VolumeMapper cropping plane or AbstractMapper clipping plane?
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

  // Create a PlaneId for image actor.
  if (this->ImageActor)
   {
   int extent[6];
   this->ImageActor->GetDisplayExtent(extent);
   this->PlaneId = 4;
   if (extent[2] == extent[3]) { this->PlaneId = 2; }
   else if (extent[0] == extent[1]) { this->PlaneId = 0; }
   }

  if (this->PlaneId >= 0)
    {
    this->GetPlaneOriginAndNormal(this->Origin, this->Normal);
    }
}

//----------------------------------------------------------------------------
void vtkPushPlaneAction::GetPlaneOriginAndNormal(double origin[3],
                                                 double normal[3])
{
  if (this->Mapper)
    {
    vtkPlane *plane =
      this->Mapper->GetClippingPlanes()->GetItem(this->PlaneId);

    plane->GetNormal(normal);
    plane->GetOrigin(origin);
    }
  else
    {
    double bounds[6];

    if (this->ImageActor)
      {
      this->ImageActor->GetDisplayBounds(bounds);
      }
    else if (this->VolumeMapper)
      {
      this->VolumeMapper->GetCroppingRegionPlanes(bounds);
      }
    else
      {
      return;
      }

    int i = this->PlaneId/2;

    normal[0] = normal[1] = normal[2] = 0.0;
    normal[i] = 1.0;

    origin[0] = 0.5*(bounds[0] + bounds[1]);
    origin[1] = 0.5*(bounds[2] + bounds[3]);
    origin[2] = 0.5*(bounds[4] + bounds[5]);
    origin[i] = bounds[this->PlaneId];

    // Transform from data coords to world coords.
    this->Transform->TransformNormal(normal, normal);
    this->Transform->TransformPoint(origin, origin);
    }
}

//----------------------------------------------------------------------------
void vtkPushPlaneAction::SetOrigin(const double o[3])
{
  double origin[3];

  this->Origin[0] = origin[0] = o[0];
  this->Origin[1] = origin[1] = o[1];
  this->Origin[2] = origin[2] = o[2];

  if (this->PlaneId < 0)
    {
    return;
    }

  if (this->Mapper)
    {
    vtkPlane *plane =
      this->Mapper->GetClippingPlanes()->GetItem(this->PlaneId);

    // Bounding checks needed!

    plane->SetOrigin(origin);
    }
  else
    {
    // Go from world coords to data coords
    this->Transform->GetInverse()->TransformPoint(origin, origin);

    int i = this->PlaneId/2;

    if (this->ImageActor)
      {
      double dataOrigin[3];
      this->ImageActor->GetInput()->GetOrigin(dataOrigin);
      double dataSpacing[3];
      this->ImageActor->GetInput()->GetSpacing(dataSpacing);
      int wholeExtent[6];
      this->ImageActor->GetInput()->GetWholeExtent(wholeExtent);
      int displayExtent[6];
      this->ImageActor->GetDisplayExtent(displayExtent);
      
      double x = (origin[i] - dataOrigin[i])/dataSpacing[i];
      if (x < wholeExtent[2*i]) { x = wholeExtent[2*i]; }
      if (x > wholeExtent[2*i+1]) { x = wholeExtent[2*i+1]; }
 
      int xi = int(floor(x));
      if ((x - xi) >= 0.5) { xi++; }

      displayExtent[2*i] = displayExtent[2*i+1] = xi; 
      this->ImageActor->SetDisplayExtent(displayExtent);
      }
    else if (this->VolumeMapper)
      {
      double region[6];
      this->VolumeMapper->GetCroppingRegionPlanes(region);
      double bounds[6];
      this->VolumeMapper->GetBounds(bounds);

      // Get the cropping plane position
      double x = origin[i];

      // Get the minimum thickness of volume allowed.
      double t = 1.0;
      vtkImageData *data = this->VolumeMapper->GetInput();
      if (data)
        {
        double spacing[3];
        data->GetSpacing(spacing);
        t = spacing[i];
        }

      // Check for collissions with the opposing plane
      if (this->PlaneId == 2*i)
        {
        if (x > region[2*i+1] - t) { x = region[2*i+1] - t; }
        if (x > bounds[2*i+1] - t) { x = bounds[2*i+1] - t; }
        }
      else
        {
        if (x < region[2*i] + t) { x = region[2*i] + t; }
        if (x < bounds[2*i] + t) { x = bounds[2*i] + t; }
        }

      // Bounding box check
      if (x < bounds[2*i]) { x = bounds[2*i]; }
      if (x > bounds[2*i+1]) { x = bounds[2*i+1]; }

      region[this->PlaneId] = x;
      this->VolumeMapper->SetCroppingRegionPlanes(region);
      }
    }
}   

//----------------------------------------------------------------------------
void vtkPushPlaneAction::SetNormal(const double normal[3])
{
  this->Normal[0] = normal[0];
  this->Normal[1] = normal[1];
  this->Normal[2] = normal[2];

  if (this->PlaneId < 0)
    {
    return;
    }

  // Do nothing for now.  Will need to implement for plane rotation.
}

