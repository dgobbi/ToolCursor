/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkSurfaceCursor.cxx,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkSurfaceCursor.h"
#include "vtkObjectFactory.h"

#include "vtkRenderer.h"
#include "vtkCamera.h"
#include "vtkActor.h"
#include "vtkVolume.h"
#include "vtkVolumeMapper.h"
#include "vtkImageActor.h"
#include "vtkProp3DCollection.h"
#include "vtkPlaneCollection.h"
#include "vtkPlane.h"
#include "vtkAssemblyPath.h"
#include "vtkProperty.h"
#include "vtkDataSetMapper.h"
#include "vtkAbstractVolumeMapper.h"
#include "vtkLookupTable.h"
#include "vtkDataSetCollection.h"
#include "vtkImageData.h"
#include "vtkPolyData.h"
#include "vtkMatrix4x4.h"
#include "vtkMath.h"
#include "vtkVolumePicker.h"
#include "vtkCommand.h"

#include "vtkSurfaceCursorShapes.h"
#include "vtkActionCursorShapes.h"
#include "vtkSurfaceCursorAction.h"
#include "vtkPushPlaneAction.h"

vtkCxxRevisionMacro(vtkSurfaceCursor, "$Revision: 1.23 $");
vtkStandardNewMacro(vtkSurfaceCursor);

//----------------------------------------------------------------------------
class vtkSurfaceCursorRenderCommand : public vtkCommand
{
public:
  static vtkSurfaceCursorRenderCommand *New(vtkSurfaceCursor *cursor) {
    return new vtkSurfaceCursorRenderCommand(cursor); };

  virtual void Execute(vtkObject *object, unsigned long, void *) {
    this->Cursor->OnRender(); };

protected:
  vtkSurfaceCursorRenderCommand(vtkSurfaceCursor *cursor) {
    this->Cursor = cursor; };

  vtkSurfaceCursor* Cursor;

private:
  static vtkSurfaceCursorRenderCommand *New(); // Not implemented.
  vtkSurfaceCursorRenderCommand(); // Not implemented.
  vtkSurfaceCursorRenderCommand(const vtkSurfaceCursorRenderCommand&);  // Not implemented.
  void operator=(const vtkSurfaceCursorRenderCommand&);  // Not implemented.
};

//----------------------------------------------------------------------------
vtkSurfaceCursor::vtkSurfaceCursor()
{
  this->DisplayPosition[0] = 0.0;
  this->DisplayPosition[1] = 0.0;

  this->Position[0] = 0.0;
  this->Position[1] = 0.0;
  this->Position[2] = 0.0;

  this->Normal[0] = 0.0;
  this->Normal[1] = 0.0;
  this->Normal[2] = 1.0;

  this->Vector[0] = 0.0;
  this->Vector[1] = 1.0;
  this->Vector[2] = 0.0;

  this->Renderer = 0;

  this->PointNormalAtCamera = 1;
  this->Modifier = 0;
  this->Mode = 0;
  this->PickFlags = 0;
  this->Shape = 0;
  this->Action = 0;
  this->ActionButton = 0;
  this->Scale = 1.0;

  this->Actor = vtkActor::New();
  this->Matrix = vtkMatrix4x4::New();
  this->Mapper = vtkDataSetMapper::New();
  this->Mapper->StaticOn();
  this->LookupTable = vtkLookupTable::New();
  this->Mapper->SetLookupTable(this->LookupTable);
  this->Mapper->UseLookupTableScalarRangeOn();
  //this->Shapes = vtkSurfaceCursorShapes::New();
  this->Actions = vtkCollection::New();
  this->Picker = vtkVolumePicker::New();

  this->LookupTable->SetRampToLinear();
  this->LookupTable->SetTableRange(0,255);
  this->LookupTable->SetNumberOfTableValues(256);
  this->LookupTable->SetSaturationRange(0,0);
  this->LookupTable->SetValueRange(0,1);
  this->LookupTable->Build();
  this->LookupTable->SetTableValue(0, 1.0, 0.0, 0.0, 1.0);
  this->LookupTable->SetTableValue(1, 0.0, 1.0, 0.0, 1.0);

  this->Actor->PickableOff();
  this->Actor->VisibilityOff();
  this->Actor->SetMapper(this->Mapper);
  this->Actor->SetUserMatrix(this->Matrix);

  vtkProperty *property = this->Actor->GetProperty();
  property->BackfaceCullingOn();

  this->Shapes = vtkActionCursorShapes::New();
  this->SetShape(0);

  // Make all the actions.  Should be done in its own subroutine.
  vtkSurfaceCursorAction *actionObject;
  actionObject = vtkSurfaceCursorAction::New();
  this->AddAction(actionObject);
  actionObject->Delete();
  actionObject = vtkPushPlaneAction::New();
  this->AddAction(actionObject);
  actionObject->Delete();

  this->RenderCommand = vtkSurfaceCursorRenderCommand::New(this);
}

//----------------------------------------------------------------------------
vtkSurfaceCursor::~vtkSurfaceCursor()
{
  this->SetRenderer(0);

  if (this->RenderCommand)
    {
    this->RenderCommand->Delete();
    }
  if (this->Matrix)
    {
    this->Matrix->Delete();
    }
  if (this->Shapes)
    {
    this->Shapes->Delete();
    }
  if (this->Actions)
    {
    this->Actions->Delete();
    }
  if (this->Mapper)
    {
    this->Mapper->Delete();
    }
  if (this->LookupTable)
    {
    this->LookupTable->Delete();
    }
  if (this->Actor)
    {
    this->Actor->Delete();
    }
  if (this->Picker)
    {
    this->Picker->Delete();
    }
}

//----------------------------------------------------------------------------
void vtkSurfaceCursor::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

//----------------------------------------------------------------------------
void vtkSurfaceCursor::SetRenderer(vtkRenderer *renderer)
{
  if (renderer == this->Renderer)
    {
    return;
    }

  if (this->Renderer)
    {
    this->Renderer->RemoveObserver(this->RenderCommand);
    this->Renderer->RemoveActor(this->Actor);
    this->Renderer->Delete();
    this->Renderer = 0;
    }

  if (renderer)
    {
    this->Renderer = renderer;
    this->Renderer->Register(this);
    this->Renderer->AddActor(this->Actor);
    this->Renderer->AddObserver(vtkCommand::StartEvent,
                                this->RenderCommand, -1);
    }

  this->Modified();
}

//----------------------------------------------------------------------------
void vtkSurfaceCursor::SetColor(int i, double r, double b, double g)
{
  if (i >= 0 && i <= 255)
    {
    double rgba[4];
    this->LookupTable->GetTableValue(i, rgba);
    if (rgba[0] != r || rgba[1] != g || rgba[2] != b)
      {
      this->LookupTable->SetTableValue(i, r, g, b, 1.0);
      this->Modified();
      }
    }
}

//----------------------------------------------------------------------------
void vtkSurfaceCursor::GetColor(int i, double rgb[3])
{
  if (i < 0) { i = 0; }
  if (i > 255) { i = 255; }

  double rgba[4];
  this->LookupTable->GetTableValue(i, rgba);

  rgb[0] = rgba[0];
  rgb[1] = rgba[1];
  rgb[2] = rgba[2];
}  

//----------------------------------------------------------------------------
void vtkSurfaceCursor::UpdatePropsForPick(vtkPicker *picker,
                                          vtkRenderer *renderer)
{
  // Go through all Prop3Ds that might be picked and update their data.
  // This is necessary if any data has changed since the last render.

  vtkPropCollection *props;
  if ( picker->GetPickFromList() )
    {
    props = picker->GetPickList();
    }
  else
    {
    props = renderer->GetViewProps();
    }

  vtkProp *prop;
  vtkCollectionSimpleIterator pit;
  props->InitTraversal(pit);
  while ( (prop = props->GetNextProp(pit)) )
    {
    vtkAssemblyPath *path;
    prop->InitPathTraversal();
    while ( (path = prop->GetNextPath()) )
      {
      if (!prop->GetPickable() || !prop->GetVisibility())
        {
        break;
        }

      vtkProp *anyProp = path->GetLastNode()->GetViewProp();
      vtkActor *actor;
      vtkVolume *volume;
      vtkImageActor *imageActor;
      
      if ( (actor = vtkActor::SafeDownCast(anyProp)) )
        {
        vtkDataSet *data = actor->GetMapper()->GetInput();
        if (data)
          {
          data->Update();
          }
        }
      else if ( (volume = vtkVolume::SafeDownCast(anyProp)) )
        {
        vtkDataSet *data = volume->GetMapper()->GetDataSetInput();
        if (data)
          {
          data->UpdateInformation();
          data->SetUpdateExtentToWholeExtent();
          data->Update();
          }
        }
      else if ( (imageActor = vtkImageActor::SafeDownCast(anyProp)) )
        {
        vtkImageData *data = imageActor->GetInput();
        if (data)
          {
          data->UpdateInformation();
          int extent[6], wextent[6], dextent[6];
          data->GetExtent(extent);
          data->GetWholeExtent(wextent);
          imageActor->GetDisplayExtent(dextent);
          if (dextent[0] == -1)
            {
            for (int i = 0; i < 6; i++) { extent[i] = wextent[i]; }
            if (extent[5] < extent[4])
              {
              extent[5] = extent[4];
              }
            }
          else
            {
            for (int i = 0; i < 3; i++)
              {
              int l = 2*i;
              int h = l+1;
              // Clip the display extent with the whole extent
              if (dextent[l] > wextent[l]) { dextent[l] = wextent[l]; }
              if (dextent[h] < wextent[h]) { dextent[h] = wextent[h]; }
              // Expand the extent to include the display extent
              if (extent[l] > dextent[l]) { extent[l] = dextent[l]; }
              if (extent[h] < dextent[h]) { extent[h] = dextent[h]; }
              }
            }
          data->SetUpdateExtent(extent);
          data->PropagateUpdateExtent();
          data->UpdateData();
          }
        }
      }
    }
}

//----------------------------------------------------------------------------
int vtkSurfaceCursor::ComputeMode(int modifier)
{
  // The PickCropping/Clipping shouldn't be hard-coded like this.

  if ( (modifier & VTK_SCURSOR_SHIFT) )
    {
    this->Picker->PickCroppingPlanesOn();
    this->Picker->PickClippingPlanesOff();
    return 1;
    }
  else if ( (modifier & VTK_SCURSOR_CONTROL) )
    {
    this->Picker->PickClippingPlanesOn();
    this->Picker->PickCroppingPlanesOff();
    return 2;
    }

  this->Picker->PickClippingPlanesOff();
  this->Picker->PickCroppingPlanesOff();
  return 0;
}

//----------------------------------------------------------------------------
int vtkSurfaceCursor::ComputeShape(int mode, int pickFlags)
{
  int shape = VTK_SCURSOR_POINTER;

  // Setting the cursor shape from the pickFlags is very crude for now
  switch (mode)
    {
    case 0:
      {
      if ((pickFlags & (VTK_SCURSOR_CLIP_PLANE | VTK_SCURSOR_CROP_PLANE)))
        {
        shape = VTK_SCURSOR_PUSHER;
        }
      else if ((pickFlags & VTK_SCURSOR_PROP3D))
        {
        shape = VTK_SCURSOR_CONE;
        }
      else
        {
        shape = VTK_SCURSOR_POINTER;
        }
      }
      break;

    case 1:
      {
      if ((pickFlags & VTK_SCURSOR_CROP_PLANE))
        {
        shape = VTK_SCURSOR_PUSHER;
        }
      else if ((pickFlags & VTK_SCURSOR_CLIP_PLANE))
        {
        shape = VTK_SCURSOR_ROCKER;
        }
      else if ((pickFlags & VTK_SCURSOR_PROP3D))
        {
        shape = VTK_SCURSOR_MOVER;
        }
      }
      break;

    case 2:
      {
      if ((pickFlags & VTK_SCURSOR_IMAGE_ACTOR) ||
          ((pickFlags & VTK_SCURSOR_VOLUME) &&
           (pickFlags & (VTK_SCURSOR_CROP_PLANE | VTK_SCURSOR_CLIP_PLANE))))
        {
        shape = VTK_SCURSOR_CROSS_SPLIT;
        }
      else if ((pickFlags & VTK_SCURSOR_VOLUME) ||
               (pickFlags & VTK_SCURSOR_ACTOR))
        {
        shape = VTK_SCURSOR_SPINNER;
        }
      }
      break;
    }

  return shape;
}

//----------------------------------------------------------------------------
int vtkSurfaceCursor::ComputeAction(int mode, int pickFlags, int button)
{
  // One button only.
  if (button != VTK_SCURSOR_B1)
    {
    return 0;
    }

  int action = 0;

  switch (mode)
    {
    case 0:
      {
      if ((pickFlags & (VTK_SCURSOR_CLIP_PLANE | VTK_SCURSOR_CROP_PLANE)))
        {
        action = VTK_SCURSOR_PUSH;
        }
      else if ((pickFlags & VTK_SCURSOR_PROP3D))
        {
        action = VTK_SCURSOR_ROTATE;
        }
      }
      break;

    case 1:
      {
      if ((pickFlags & VTK_SCURSOR_CROP_PLANE))
        {
        action = VTK_SCURSOR_PUSH;
        }
      else if ((pickFlags & VTK_SCURSOR_PROP3D))
        {
        action = VTK_SCURSOR_SPIN;
        }
      }
      break;
    }

  return action;
}

//----------------------------------------------------------------------------
int vtkSurfaceCursor::ComputePickFlags(vtkVolumePicker *picker)
{
  const double planeTol = 1e-6;
  const double normalTol = 1e-15;

  int pickFlags = 0;

  vtkProp3DCollection *props = picker->GetProp3Ds();
  vtkCollectionSimpleIterator pit;
  props->InitTraversal(pit);
  vtkProp3D *prop = props->GetNextProp3D(pit);
  vtkAbstractMapper3D *mapper = picker->GetMapper();

  if (!prop)
    {
    // No prop, nothing to do
    return 0;
    }

  if (mapper && picker->GetClippingPlaneId() >= 0)
    {
    // Make sure that our Position lies on the plane and that our Normal
    // is perpendicular to the plane.

    vtkPlane *plane = mapper->GetClippingPlanes()->GetItem(
      picker->GetClippingPlaneId());

    double u[3];
    vtkMath::Cross(plane->GetNormal(), picker->GetPickNormal(), u);

    if (fabs(plane->EvaluateFunction(picker->GetPickPosition())) < planeTol &&
        vtkMath::Norm(u) < normalTol)
      {
      pickFlags = (pickFlags | VTK_SCURSOR_CLIP_PLANE);
      }
    }

  if (mapper->IsA("vtkVolumeMapper") && picker->GetCroppingPlaneId() >= 0)
    {
    // Also ensure that our Position lies on the cropping plane
    int planeId = picker->GetCroppingPlaneId();
    vtkVolumeMapper *volumeMapper = static_cast<vtkVolumeMapper *>(mapper); 

    double bounds[6];
    volumeMapper->GetCroppingRegionPlanes(bounds);

    double mapperPos[3];
    picker->GetMapperPosition(mapperPos);

    double planeNormal[3], mapperNormal[3], u[3];
    planeNormal[0] = planeNormal[1] = planeNormal[2] = 0.0;
    planeNormal[planeId/2] = 1.0;
    picker->GetMapperNormal(mapperNormal);
    vtkMath::Cross(planeNormal, mapperNormal, u);

    if (fabs(mapperPos[planeId/2] - bounds[planeId]) < planeTol &&
        vtkMath::Norm(u) < normalTol)
      {
      pickFlags = (pickFlags | VTK_SCURSOR_CROP_PLANE);
      }
    }

  if (prop->IsA("vtkActor"))
    {
    pickFlags = (pickFlags | VTK_SCURSOR_ACTOR);
    }
  else if (prop->IsA("vtkVolume"))
    {
    pickFlags = (pickFlags | VTK_SCURSOR_VOLUME);
    }
  else if (prop->IsA("vtkImageActor"))
    {
    pickFlags = (pickFlags | VTK_SCURSOR_IMAGE_ACTOR);
    }

  return pickFlags;
}

//----------------------------------------------------------------------------
void vtkSurfaceCursor::ComputePosition()
{
  if (!this->Renderer)
    {
    return;
    }

  int x = this->DisplayPosition[0];
  int y = this->DisplayPosition[1];

  // Update the props that might be picked.  This is necessary
  // if there hasn't been a Render since the last change.
  this->UpdatePropsForPick(this->Picker, this->Renderer);

  // Do the pick!
  vtkVolumePicker *picker = this->Picker;
  picker->Pick(x, y, 0, this->Renderer);
  picker->GetPickPosition(this->Position);
  picker->GetPickNormal(this->Normal);

  // Direct the normal towards the camera if PointNormalAtCamera is On
  if (this->PointNormalAtCamera &&
      vtkMath::Dot(this->Renderer->GetActiveCamera()
                   ->GetDirectionOfProjection(),
                   this->Normal) > 0)
    {
    this->Normal[0] = -this->Normal[0];
    this->Normal[1] = -this->Normal[1];
    this->Normal[2] = -this->Normal[2];
    }

  // Compute the "state" from the picked information.  The "state" needs to
  // be split into two ints: there should be a "pick info" bitfield that
  // describes what's under the cursor, and another "state" that describes the
  // ongoing interaction and which is locked until the action is completed. 
  this->PickFlags = this->ComputePickFlags(this->Picker);

  // Compute the cursor shape from the state.  
  this->SetShape(this->ComputeShape(this->Mode, this->PickFlags));

  // Compute an "up" vector for the cursor.
  this->ComputeVectorFromNormal(this->Normal, this->Vector,
                                this->Mapper, this->Renderer,
                                this->Shapes->GetShapeFlags(this->Shape));

  // Compute the pose matrix for the cursor
  this->ComputeMatrix(this->Position, this->Normal, this->Vector,
                      this->Matrix);

  // Scale for the cursor to always be the same number of pixels across.
  double scale = this->ComputeScale(this->Position, this->Renderer);

  this->Actor->SetScale(scale*this->Scale);
}

//----------------------------------------------------------------------------
void vtkSurfaceCursor::Modified()
{
  this->Mapper->Modified();
  this->Superclass::Modified();
}

//----------------------------------------------------------------------------
void vtkSurfaceCursor::OnRender()
{
  // Compute the position when the Renderer renders, since it needs to
  // update all of the props in the scene.
  this->ComputePosition();
  // Don't show cursor if nothing is underneath of it.
  int visibility = (this->MouseInRenderer != 0 && this->PickFlags != 0);
  this->Actor->SetVisibility(visibility);
}

//----------------------------------------------------------------------------
void vtkSurfaceCursor::MoveToDisplayPosition(double x, double y)
{
  this->SetDisplayPosition(x, y);

  if (this->Action)
    {
    vtkSurfaceCursorAction *actionObject =
      static_cast<vtkSurfaceCursorAction *>(
        this->Actions->GetItemAsObject(this->Action));

    if (actionObject)
      {
      actionObject->DoAction();
      }
   }
}

//----------------------------------------------------------------------------
int vtkSurfaceCursor::GetRequestingFocus()
{
  // Right now the cursor is hard-coded to take focus for clipping planes
  return ((this->PickFlags & VTK_SCURSOR_CLIP_PLANE) ||
          (this->PickFlags & VTK_SCURSOR_CROP_PLANE));
}

//----------------------------------------------------------------------------
int vtkSurfaceCursor::GetVisibility()
{
  return this->Actor->GetVisibility();
}

//----------------------------------------------------------------------------
void vtkSurfaceCursor::SetMouseInRenderer(int inside)
{
  if (this->MouseInRenderer == inside)
    {
    return;
    }
    
  this->MouseInRenderer = inside;
  this->Modified();
}

//----------------------------------------------------------------------------
void vtkSurfaceCursor::SetModifier(int modifier)
{
  if (this->Modifier == modifier)
    {
    return;
    }
    
  // Map the modifier to the Mode.
  this->SetMode(this->ComputeMode(modifier));

  // Map the modifier to the Action.
  // Do an XOR to find out what bits have changed.
  int bitsChanged = (this->Modifier ^ modifier);
  // Check for mouse button changes.
  if (bitsChanged & (VTK_SCURSOR_B1 | VTK_SCURSOR_B2 | VTK_SCURSOR_B3))
    {
    int bitsSet = (bitsChanged & modifier);
    int bitsCleared = (bitsChanged & ~modifier);

    // Check whether the active action button was released.
    if (this->Action && (bitsCleared & this->ActionButton))
      {
      this->SetAction(0);
      this->ActionButton = 0;
      }
    else if (!this->Action && bitsSet)
      {
      int button = 0;
      if ((bitsSet & VTK_SCURSOR_B1)) { button = VTK_SCURSOR_B1; }
      else if ((bitsSet & VTK_SCURSOR_B2)) { button = VTK_SCURSOR_B2; }
      else if ((bitsSet & VTK_SCURSOR_B3)) { button = VTK_SCURSOR_B3; }
      this->SetAction(this->ComputeAction(this->Mode, this->PickFlags,
                                          button));
      this->ActionButton = button;
      }
    }

  this->Modifier = modifier;
  this->Modified();
}

//----------------------------------------------------------------------------
void vtkSurfaceCursor::SetAction(int action)
{
  if (action == this->Action)
    {
    return;
    }

  if (this->Action)
    {
    vtkSurfaceCursorAction *actionObject =
      static_cast<vtkSurfaceCursorAction *>(
        this->Actions->GetItemAsObject(this->Action));

    if (actionObject)
      {
      actionObject->StopAction();
      actionObject->SetSurfaceCursor(0);
      }
    }

  this->Action = action;
  this->Modified();

  if (action)
    {
    vtkSurfaceCursorAction *actionObject =
      static_cast<vtkSurfaceCursorAction *>(
        this->Actions->GetItemAsObject(this->Action));

    if (actionObject)
      {
      actionObject->SetSurfaceCursor(this);
      actionObject->StartAction();
      }
    }
}

//----------------------------------------------------------------------------
int vtkSurfaceCursor::AddAction(vtkSurfaceCursorAction *action)
{
  this->Actions->AddItem(action);
  
  return (this->Actions->GetNumberOfItems() - 1);
}
 
//----------------------------------------------------------------------------
void vtkSurfaceCursor::SetMode(int mode)
{
  if (this->Mode == mode)
    {
    return;
    }
    
  this->Mode = mode;
  this->Modified();
}

//----------------------------------------------------------------------------
void vtkSurfaceCursor::SetShape(int shape)
{
  vtkDataSet *data = this->Shapes->GetShapeData(shape);

  if (data)
    {
    this->Mapper->SetInput(data);
    this->Shape = shape;
    }
}

//----------------------------------------------------------------------------
double vtkSurfaceCursor::ComputeScale(const double position[3],
                                      vtkRenderer *renderer)
{
  // Find the cursor scale factor such that 1 data unit length
  // equals 1 screen pixel at the cursor's distance from the camera.
  // Start by computing the height of the window at the cursor position.
  double worldHeight = 1.0;
  vtkCamera *camera = renderer->GetActiveCamera();
  if (camera->GetParallelProjection())
    {
    worldHeight = 2*camera->GetParallelScale();
    }
  else
    {
    double cameraPosition[3];
    camera->GetPosition(cameraPosition);
    worldHeight = 2*(sqrt(vtkMath::Distance2BetweenPoints(position,
                                                        cameraPosition))
                     * tan(0.5*camera->GetViewAngle()/57.296));
    }

  // Compare world height to window height.
  int windowHeight = renderer->GetSize()[1];
  double scale = 1.0;
  if (windowHeight > 0)
    {
    scale = worldHeight/windowHeight;
    }

  return scale;
}

//----------------------------------------------------------------------------
void vtkSurfaceCursor::ComputeMatrix(const double p[3], const double n[3],
                                     const double v[3], vtkMatrix4x4 *matrix)
{
  double u[3];
  vtkMath::Cross(v, n, u);

  for (int j = 0; j < 3; j++)
    {
    matrix->SetElement(j, 0, u[j]);
    matrix->SetElement(j, 1, v[j]);
    matrix->SetElement(j, 2, n[j]);
    matrix->SetElement(j, 3, p[j]);
    }
  matrix->Modified();
}

//----------------------------------------------------------------------------
void vtkSurfaceCursor::ComputeVectorFromNormal(const double normal[3],
                                               double vector[3],
                                               vtkDataSetMapper *cursorMapper,
                                               vtkRenderer *renderer,
                                               int cursorFlags)
{
  // Get the camera orientation
  vtkCamera *camera = renderer->GetActiveCamera();
  vtkMatrix4x4 *matrix = camera->GetViewTransformMatrix();

  // These ints say how we want to create the vector
  int direction = 1; // We want to align the cursor y vector with...
  int primary = 1;   // the camera view up vector if possible...
  int secondary = 2; // or with the view plane normal otherwise.

  // If the data is "flat" and the flat dimension is not the z dimension,
  // then point the flat side at the camera.
  if ( (cursorFlags & VTK_SCURSOR_FLATX) )
    {
    direction = 0;
    primary = 2;
    secondary = 1;
    }
  else if ( (cursorFlags & VTK_SCURSOR_FLATY) )
    {
    direction = 1;
    primary = 2;
    secondary = 1;
    }

  // Get primary direction from camera, orthogonalize to the normal.
  double u[3];
  u[0] = matrix->GetElement(primary,0);
  u[1] = matrix->GetElement(primary,1);
  u[2] = matrix->GetElement(primary,2);

  // dot will be 1.0 if primary and normal are the same
  double dot = vtkMath::Dot(normal, u);

  if (dot > 0.999)
    {
    // blend the vector with the secondary for stability
    u[0] += matrix->GetElement(secondary,0) * (dot - 0.999);
    u[1] += matrix->GetElement(secondary,1) * (dot - 0.999);
    u[2] += matrix->GetElement(secondary,2) * (dot - 0.999);
    }

  vtkMath::Cross(normal, u, u);
  if (direction == 1)
    {
    vtkMath::Cross(u, normal, u); 
    }

  double norm = vtkMath::Norm(u);
  vector[0] = u[0]/norm;
  vector[1] = u[1]/norm;
  vector[2] = u[2]/norm;
}

