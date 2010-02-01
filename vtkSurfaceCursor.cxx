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
#include "vtkGeometricCursorShapes.h"
#include "vtkSurfaceCursorAction.h"
#include "vtkPushPlaneAction.h"
#include "vtkPanCameraAction.h"
#include "vtkRotateCameraAction.h"
#include "vtkSpinCameraAction.h"
#include "vtkZoomCameraAction.h"

vtkCxxRevisionMacro(vtkSurfaceCursor, "$Revision: 1.44 $");
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
  this->ActionButtons = 0;
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
  this->Shapes = vtkSurfaceCursorShapes::New();
  this->Actions = vtkCollection::New();
  this->ShapeBindings = vtkIntArray::New();
  this->ShapeBindings->SetName("ShapeBindings");
  this->ShapeBindings->SetNumberOfComponents(4);
  this->ActionBindings = vtkIntArray::New();
  this->ActionBindings->SetName("ActionBindings");
  this->ActionBindings->SetNumberOfComponents(4);
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

  // Insert a null action, so that actions start at 1
  vtkSurfaceCursorAction *action = vtkSurfaceCursorAction::New();
  this->Actions->AddItem(action);
  action->Delete();

  // Insert a blank shape, so that shapes start at 1
  // (shape 0 will display the system cursor)
  vtkDataSet *data = vtkPolyData::New();
  this->Shapes->AddShape("", data, 0);
  data->Delete();

  this->RenderCommand = vtkSurfaceCursorRenderCommand::New(this);
}

//----------------------------------------------------------------------------
vtkSurfaceCursor::~vtkSurfaceCursor()
{
  this->SetRenderer(0);

  if (this->RenderCommand) { this->RenderCommand->Delete(); }
  if (this->Matrix) { this->Matrix->Delete(); }
  if (this->Shapes) { this->Shapes->Delete(); }
  if (this->Actions) { this->Actions->Delete(); }
  if (this->ShapeBindings) { this->ShapeBindings->Delete(); }
  if (this->ActionBindings) { this->ActionBindings->Delete(); }
  if (this->Mapper) { this->Mapper->Delete(); }
  if (this->LookupTable) { this->LookupTable->Delete(); }
  if (this->Actor) { this->Actor->Delete(); }
  if (this->Picker) { this->Picker->Delete(); }
}

//----------------------------------------------------------------------------
void vtkSurfaceCursor::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

//----------------------------------------------------------------------------
void vtkSurfaceCursor::BindDefaultActions()
{
  vtkSurfaceCursorShapes *actionShapes = vtkActionCursorShapes::New();
  vtkSurfaceCursorShapes *geometricShapes = vtkGeometricCursorShapes::New();
  vtkSurfaceCursorAction *pushAction = vtkPushPlaneAction::New();
  vtkSurfaceCursorAction *rotateAction = vtkRotateCameraAction::New();
  vtkSurfaceCursorAction *spinAction = vtkSpinCameraAction::New();
  vtkSurfaceCursorAction *panAction = vtkPanCameraAction::New();
  vtkSurfaceCursorAction *zoomAction = vtkZoomCameraAction::New();

  int action, shape, mode, modifier, pickInfo;

  // ============ All Bindings for Mode 0 =============
  mode = 0;

  // In terms of the "PickInfo" flags, it is necessary to start at the
  // default bindings and work up to more specific collections of
  // PickInfo flags.

  // ------------ Default Bindings ------------- 
  pickInfo = 0;

  // Binding for "Spinner" cursor and action, when user clicks background
  modifier = 0;
  shape = this->AddShape(actionShapes, "Spin");
  action = this->AddAction(spinAction);
  this->BindShape(shape, mode, pickInfo, modifier | VTK_SCURSOR_B1);
  this->BindAction(action, mode, pickInfo, modifier | VTK_SCURSOR_B1);

  // Binding for "Mover" cursor and action
  modifier = VTK_SCURSOR_SHIFT;
  shape = this->AddShape(actionShapes, "Move");
  action = this->AddAction(panAction);
  this->BindShape(shape, mode, pickInfo, modifier);
  this->BindShape(shape, mode, pickInfo, modifier | VTK_SCURSOR_B1);
  this->BindAction(action, mode, pickInfo, modifier | VTK_SCURSOR_B1);

  // Also bind to middle button
  modifier = 0;
  this->BindShape(shape, mode, pickInfo, modifier | VTK_SCURSOR_B3);
  this->BindAction(action, mode, pickInfo, modifier | VTK_SCURSOR_B3);

  // Binding for "Zoom" cursor and action
  modifier = VTK_SCURSOR_CONTROL;
  shape = this->AddShape(actionShapes, "Zoom");
  action = this->AddAction(zoomAction);
  this->BindShape(shape, mode, pickInfo, modifier);
  this->BindShape(shape, mode, pickInfo, modifier | VTK_SCURSOR_B1);
  this->BindAction(action, mode, pickInfo, modifier | VTK_SCURSOR_B1);

  // Also bind to right button
  modifier = 0;
  this->BindShape(shape, mode, pickInfo, modifier | VTK_SCURSOR_B2);
  this->BindAction(action, mode, pickInfo, modifier | VTK_SCURSOR_B2);

  // ------------ Default Prop3D Bindings------------- 
  pickInfo = VTK_SCURSOR_PROP3D;

  // Binding for "Cone" cursor
  modifier = 0;
  shape = this->AddShape(geometricShapes, "Cone");
  this->BindShape(shape, mode, pickInfo, modifier);

  // Binding for "Rotator" cursor and action
  modifier = 0;
  shape = this->AddShape(actionShapes, "Rotate");
  action = this->AddAction(rotateAction);
  this->BindShape(shape, mode, pickInfo, modifier | VTK_SCURSOR_B1);
  this->BindAction(action, mode, pickInfo, modifier | VTK_SCURSOR_B1);

  // ------------ Bindings for Prop3Ds with Clip/Crop Planes -------------
  pickInfo = (VTK_SCURSOR_PROP3D |
              VTK_SCURSOR_CLIP_PLANE |
              VTK_SCURSOR_CROP_PLANE);

  // Binding for "SplitCross" cursor
  modifier = 0;
  shape = this->AddShape(geometricShapes, "SplitCross");
  this->BindShape(shape, mode, pickInfo, modifier);

  // Binding for "PushPlane" with "Pusher" cursor
  modifier = 0;
  shape = this->AddShape(actionShapes, "Push");
  action = this->AddAction(pushAction);
  this->BindShape(shape, mode, pickInfo, modifier | VTK_SCURSOR_B1);
  this->BindAction(action, mode, pickInfo, modifier | VTK_SCURSOR_B1);


  actionShapes->Delete();
  geometricShapes->Delete();
  rotateAction->Delete();
  pushAction->Delete();
  spinAction->Delete();
  panAction->Delete();
  zoomAction->Delete();
}

//----------------------------------------------------------------------------
void vtkSurfaceCursor::BindShape(int shape, int mode,
                                 int pickFlags, int modifier)
{
  this->AddBinding(this->ShapeBindings, shape, mode, pickFlags, modifier);
}

//----------------------------------------------------------------------------
void vtkSurfaceCursor::BindAction(int action, int mode,
                                  int pickFlags, int modifier)
{
  this->AddBinding(this->ActionBindings, action, mode, pickFlags, modifier);
}

//----------------------------------------------------------------------------
int vtkSurfaceCursor::FindShape(int mode, int pickFlags, int modifier)
{
  // Search the shape bindings and return the first match

  int i = this->ResolveBinding(this->ShapeBindings, 0,
                               mode, pickFlags, modifier);
  if (i < 0)
    {
    return 0;
    }

  int tuple[4];
  this->ShapeBindings->GetTupleValue(i, tuple);

  return tuple[3];
}

//----------------------------------------------------------------------------
int vtkSurfaceCursor::FindAction(int mode, int pickFlags, int modifier)
{
  // Search the action bindings and return the first match

  int i = this->ResolveBinding(this->ActionBindings, 0,
                               mode, pickFlags, modifier);
  if (i < 0)
    {
    return 0;
    }

  int tuple[4];
  this->ActionBindings->GetTupleValue(i, tuple);

  return tuple[3];
}

//----------------------------------------------------------------------------
int vtkSurfaceCursor::FindActionButtons(int mode, int pickFlags, int modifier)
{
  // Find all matching actions and return a bitmask of the mouse buttons
  // for those actions.

  int modifierMask = 0;
  modifier |= VTK_SCURSOR_BUTTON_MASK;

  int tuple[4];
  int j = 0;
  for (;;)
    {
    int i = this->ResolveBinding(this->ActionBindings, j,
                                 mode, pickFlags, modifier);
    if (i < 0) { break; }

    this->ActionBindings->GetTupleValue(i, tuple);
    modifierMask |= tuple[2];

    j = i + 1;
    }

  return (modifierMask & VTK_SCURSOR_BUTTON_MASK);
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
int vtkSurfaceCursor::ComputePickFlags(vtkVolumePicker *picker)
{
  const double planeTol = 1e-6;
  const double normalTol = 1e-15;

  int pickFlags = 0;

  vtkProp3D *prop = picker->GetProp3D();
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

  if (mapper && mapper->IsA("vtkVolumeMapper") &&
      picker->GetCroppingPlaneId() >= 0)
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

  if (mapper && mapper->IsA("vtkMapper"))
    {
    pickFlags = (pickFlags | VTK_SCURSOR_ACTOR);
    }
  else if (mapper && mapper->IsA("vtkAbstractVolumeMapper"))
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

  // Allow the action to constrain the cursor
  if (this->Action)
    {
    vtkSurfaceCursorAction *actionObject =
      static_cast<vtkSurfaceCursorAction *>(
        this->Actions->GetItemAsObject(this->Action));

    if (actionObject)
      {
      actionObject->ConstrainCursor(this->Position, this->Normal);
      }
    }

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

  // Check to see if the PickFlags have changed
  int pickFlags = this->ComputePickFlags(this->Picker);
  if (pickFlags != this->PickFlags &&
      (this->Modifier & VTK_SCURSOR_BUTTON_MASK) == 0 &&
      !this->Action)
    {
    // Compute the cursor shape from the state.  
    this->SetShape(this->FindShape(this->Mode, pickFlags, this->Modifier));

    // See if we need focus for potential button actions
    this->ActionButtons = this->FindActionButtons(this->Mode, pickFlags,
                                                 this->Modifier);
    }

  this->PickFlags = pickFlags;

  // Compute an "up" vector for the cursor.
  this->ComputeVectorFromNormal(this->Position, this->Normal, this->Vector,
                                this->Renderer,
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
  int visibility = (this->IsInViewport != 0 && this->Shape != 0);
  this->Actor->SetVisibility(visibility);
}

//----------------------------------------------------------------------------
void vtkSurfaceCursor::SetDisplayPosition(double x, double y)
{
  if (this->DisplayPosition[0] == x && this->DisplayPosition[1] == y)
    {
    return;
    }

  this->DisplayPosition[0] = x;
  this->DisplayPosition[1] = y;

  this->Modified();
}

//----------------------------------------------------------------------------
void vtkSurfaceCursor::MoveToDisplayPosition(double x, double y)
{
  this->SetDisplayPosition(x, y);

  if (this->Renderer)
    {
    this->SetIsInViewport(this->Renderer->IsInViewport(x, y));
    }

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
int vtkSurfaceCursor::GetVisibility()
{
  return this->Actor->GetVisibility();
}

//----------------------------------------------------------------------------
void vtkSurfaceCursor::SetIsInViewport(int inside)
{
  if (this->IsInViewport == inside)
    {
    return;
    }
    
  this->IsInViewport = inside;
  this->Modified();
}

//----------------------------------------------------------------------------
int vtkSurfaceCursor::ButtonBit(int button)
{
  static int buttonmap[6] = { 0, VTK_SCURSOR_B1, VTK_SCURSOR_B2,
    VTK_SCURSOR_B3, VTK_SCURSOR_B4, VTK_SCURSOR_B5 };

  if (button <= 0 || button > 5)
    {
    return 0;
    }

  return buttonmap[button];
}

//----------------------------------------------------------------------------
int vtkSurfaceCursor::PressButton(int button)
{
  int buttonBit = this->ButtonBit(button);
  int started = 0;

  this->SetModifier(this->Modifier | buttonBit);

  if (buttonBit && !this->Action)
    {
    this->SetAction(this->FindAction(this->Mode, this->PickFlags,
                                     this->Modifier));
    if (this->Action)
      {
      this->ActionButton = button;
      started = 1;
      }
    }

  return started;
}

//----------------------------------------------------------------------------
int vtkSurfaceCursor::ReleaseButton(int button)
{
  int buttonBit = this->ButtonBit(button);
  int concluded = 0;

  if (this->Action && button == this->ActionButton)
    {
    this->SetAction(0);
    this->ActionButton = 0;
    concluded = 1;
    }

  // Force SetModifier to see a change
  this->Modifier = (this->Modifier | buttonBit);
  this->SetModifier(this->Modifier & ~buttonBit);

  return concluded;
}

//----------------------------------------------------------------------------
void vtkSurfaceCursor::SetModifier(int modifier)
{
  if (this->Modifier == modifier)
    {
    return;
    }

  if ( ((this->Modifier & VTK_SCURSOR_BUTTON_MASK) == 0 ||
        (modifier & VTK_SCURSOR_BUTTON_MASK) == 0) &&
       !this->Action)
    {
    // Compute the cursor shape from the state.  
    this->SetShape(this->FindShape(this->Mode, this->PickFlags, modifier));

    // See if we need focus for potential button actions
    this->ActionButtons = this->FindActionButtons(this->Mode, this->PickFlags,
                                                  modifier);
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
int vtkSurfaceCursor::AddShape(vtkSurfaceCursorShapes *shapes,
                               const char *name)
{
  int i = shapes->GetShapeIndex(name);
  if (i < 0)
    {
    vtkErrorMacro("The specified shape \"" << name << "\" is not in "
                   << shapes->GetClassName());
    return -1;
    }

  this->Shapes->AddShape(name, shapes->GetShapeData(i),
                         shapes->GetShapeFlags(i));
  
  return (this->Shapes->GetNumberOfShapes() - 1);
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
    vtkMatrix4x4 *matrix = camera->GetViewTransformMatrix();
    // Get a 3x3 matrix with the camera orientation
    double cvz[3];
    cvz[0] = matrix->GetElement(2, 0);
    cvz[1] = matrix->GetElement(2, 1);
    cvz[2] = matrix->GetElement(2, 2);

    double cameraPosition[3];
    camera->GetPosition(cameraPosition);

    double v[3];
    v[0] = cameraPosition[0] - position[0];
    v[1] = cameraPosition[1] - position[1];
    v[2] = cameraPosition[2] - position[2];

    worldHeight = 2*(vtkMath::Dot(v,cvz)
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
void vtkSurfaceCursor::ComputeVectorFromNormal(const double position[3],
                                               const double normal[3],
                                               double vector[3],
                                               vtkRenderer *renderer,
                                               int cursorFlags)
{
  // Mask out all but the relevant flags
  cursorFlags = (cursorFlags & VTK_SCURSOR_ORIENT);

  // Get the camera orientation
  vtkCamera *camera = renderer->GetActiveCamera();
  vtkMatrix4x4 *matrix = camera->GetViewTransformMatrix();

  // Get a 3x3 matrix with the camera orientation
  double mat[3][3];
  for (int i = 0; i < 3; i++)
    {
    mat[i][0] = matrix->GetElement(i, 0);
    mat[i][1] = matrix->GetElement(i, 1);
    mat[i][2] = matrix->GetElement(i, 2);
    }

  if (cursorFlags == VTK_SCURSOR_RADIALX ||
      cursorFlags == VTK_SCURSOR_RADIALY)
    {
    // Make "y" point away from camera axis
    double f[3];
    camera->GetFocalPoint(f);

    double v[3];
    v[0] = position[0] - f[0];
    v[1] = position[1] - f[1];
    v[2] = position[2] - f[2];

    // Make "x" perpendicular to new "y"
    vtkMath::Cross(v, mat[2], mat[0]);
    vtkMath::Normalize(mat[0]);

    // Orthogonalize "y" wrt to camera axis
    vtkMath::Cross(mat[2], mat[0], mat[1]);
    vtkMath::Normalize(mat[1]);
    }

  // These ints say how we want to create the vector
  int direction = 1; // We want to align the cursor y vector with...
  int primary = 1;   // the camera's y vector if possible...
  int secondary = 2; // or with the camera's z vector otherwise.

  // If the data is "flat" and the flat dimension is not the z dimension,
  // then point the flat side at the camera.
  if (cursorFlags == VTK_SCURSOR_FLATX)
    {
    direction = 0;
    primary = 2;
    secondary = 1;
    }
  else if (cursorFlags == VTK_SCURSOR_FLATY)
    {
    direction = 1;
    primary = 2;
    secondary = 1;
    }
  else if (cursorFlags == VTK_SCURSOR_RADIALX)
    {
    direction = 0;
    primary = 1;
    secondary = 0;
    }
  else if (cursorFlags == VTK_SCURSOR_RADIALY)
    {
    direction = 1;
    primary = 1;
    secondary = 0;
    }

  // Get primary direction from camera, orthogonalize to the normal.
  double u[3];
  u[0] = mat[primary][0];
  u[1] = mat[primary][1];
  u[2] = mat[primary][2];

  // dot will be 1.0 if primary and normal are the same
  double dot = vtkMath::Dot(normal, u);

  if (dot > 0.999)
    {
    // blend the vector with the secondary for stability
    u[0] += mat[secondary][0] * (dot - 0.999);
    u[1] += mat[secondary][1] * (dot - 0.999);
    u[2] += mat[secondary][2] * (dot - 0.999);
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

//----------------------------------------------------------------------------
// A couple useful diagnostic routines
/*
static void PrintModifier(int modifier);
static void PrintFlags(int flags);
*/

//----------------------------------------------------------------------------
void vtkSurfaceCursor::AddBinding(vtkIntArray *array, int item, int mode,
                                  int pickFlags, int modifier)
{
  // The bindings are sorted by mode and by the number of modifier bits.
  // Exact matches replace the previous entry.

  //cerr << "AddBinding " << array->GetName() << " ";
  //PrintFlags(pickFlags);
  //cerr << " ";
  //PrintModifier(modifier);
  //cerr << " " << item << "\n";

  int i = vtkSurfaceCursor::ResolveBinding(array, 0,
                                           mode, pickFlags, modifier);
  int n = array->GetNumberOfTuples();

  int tuple[4];
  tuple[0] = 0;
  tuple[1] = 0;
  tuple[2] = 0;
  tuple[3] = 0;

  if (i >= 0)
    {
    // If it's an exact match, then replace
    array->GetTupleValue(i, tuple);
    if (tuple[0] == mode && tuple[1] == pickFlags && tuple[2] == modifier)
      {
      tuple[3] = item;
      array->SetTupleValue(i, tuple);
      return;
      }
    }
  else
    {
    // If it wasn't resolved, then we append to the end of the list
    i = n;
    }

  //cerr << "inserting ";
  //cerr << item << " (" << i << " of " << (n+1) << "): " << mode << " ";
  //PrintFlags(pickFlags);
  //cerr << " ";
  //PrintModifier(modifier);
  //cerr << "\n";

  // Extend array by one value.  Actually, this is just a dummy tuple,
  // it's just the easiest way of increasing the size of the array.
  array->InsertNextTupleValue(tuple);

  //cerr << "n = " << n << " i = " << i << "\n";
  // Shuffle values up by one
  for (int j = n; j > i; j--)
    {
    array->GetTupleValue(j-1, tuple);
    array->SetTupleValue(j, tuple);
    }
 
  // Set the tuple at the desired index
  tuple[0] = mode;
  tuple[1] = pickFlags;
  tuple[2] = modifier;
  tuple[3] = item;

  array->SetTupleValue(i, tuple);
}

//----------------------------------------------------------------------------
int vtkSurfaceCursor::ResolveBinding(vtkIntArray *array, int start,
                                     int mode, int pickFlags, int modifier)
{
  // The following rules are used to resolve a binding:
  //
  // 1) the mode must match the binding exactly
  //
  // 2) all bits of each of the PROP3D and PLANE sections of the pick
  //    flags must have a matching bit in the binding, unless that section
  //    of the binding is empty.
  //
  // 3) the modifier bits must include all bits in the binding.

  //cerr << "ResolveBinding " << array->GetName() << " ";
  //PrintFlags(pickFlags);
  //cerr << " ";
  //PrintModifier(modifier);
  //cerr << "\n";

  int propType = (pickFlags & VTK_SCURSOR_PROP3D);
  int planeType = (pickFlags & VTK_SCURSOR_PLANE);

  int tuple[4];
  int n = array->GetNumberOfTuples();
  for (int i = start; i < n; i++)
    {
    array->GetTupleValue(i, tuple);

    //cerr << "item " << i << " of " << n << ": " << tuple[0] << " ";
    //PrintFlags(tuple[1]);
    //cerr << " ";
    //PrintModifier(tuple[2]);
    //cerr << "\n";

    if (tuple[0] == mode)
      {
      if (((tuple[1] & VTK_SCURSOR_PROP3D) == 0 ||
           (propType != 0 && (tuple[1] & propType) == propType)) &&
          ((tuple[1] & VTK_SCURSOR_PLANE) == 0 ||
           (planeType != 0 && (tuple[1] & planeType) == planeType)))
        {
        if ((tuple[2] & modifier) == tuple[2])
          {
          //cerr << "resolved " << tuple[3] << "\n";
          return i;
          }
        }
      }
    else if (tuple[0] > mode)
      {
      break;
      }
    }

  //cerr << "not resolved\n";
  return -1;
}   

/*
static void PrintModifier(int modifier)
{
  cerr << "( ";
  if (modifier & VTK_SCURSOR_SHIFT) { cerr << "SHIFT "; }
  if (modifier & VTK_SCURSOR_CAPS) { cerr << "CAPS "; }
  if (modifier & VTK_SCURSOR_CONTROL) { cerr << "CONTROL "; }
  if (modifier & VTK_SCURSOR_META) { cerr << "META "; }
  if (modifier & VTK_SCURSOR_ALT) { cerr << "ALT "; }
  if (modifier & VTK_SCURSOR_B1) { cerr << "B1 "; }
  if (modifier & VTK_SCURSOR_B2) { cerr << "B2 "; }
  if (modifier & VTK_SCURSOR_B3) { cerr << "B3 "; }
  cerr << ")";
}

static void PrintFlags(int flags)
{
  cerr << "( ";
  if ((flags & VTK_SCURSOR_PROP3D) == VTK_SCURSOR_PROP3D)
    {
    cerr << "PROP3DS ";
    }
  else
    {
    if (flags & VTK_SCURSOR_ACTOR) { cerr << "ACTOR "; }
    if (flags & VTK_SCURSOR_VOLUME) { cerr << "VOLUME "; }
    if (flags & VTK_SCURSOR_IMAGE_ACTOR) { cerr << "IMAGE_ACTOR "; }
    }
  if ((flags & VTK_SCURSOR_CLIP_PLANE) && (flags & VTK_SCURSOR_CROP_PLANE))
    {
    cerr << "PLANES ";
    }
  else
    {
    if (flags & VTK_SCURSOR_CLIP_PLANE) { cerr << "CLIP_PLANE "; }
    if (flags & VTK_SCURSOR_CROP_PLANE) { cerr << "CROP_PLANE "; }
    }
  cerr << ")";
}
*/

