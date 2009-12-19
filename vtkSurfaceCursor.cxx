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
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
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
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkDataArray.h"
#include "vtkDoubleArray.h"
#include "vtkIntArray.h"
#include "vtkPointData.h"
#include "vtkMatrix4x4.h"
#include "vtkMath.h"
#include "vtkVolumePicker.h"
#include "vtkCommand.h"

#include "vtkImplicitModeller.h"
#include "vtkContourFilter.h"
#include "vtkStripper.h"
#include "vtkPolyDataNormals.h"
#include "vtkReverseSense.h"
#include "vtkWarpTo.h"
#include "vtkTransform.h"

vtkCxxRevisionMacro(vtkSurfaceCursor, "$Revision: 1.16 $");
vtkStandardNewMacro(vtkSurfaceCursor);

//----------------------------------------------------------------------------
class vtkSurfaceCursorCommand : public vtkCommand
{
public:
  static vtkSurfaceCursorCommand *New(vtkSurfaceCursor *cursor) {
    return new vtkSurfaceCursorCommand(cursor); };

  virtual void Execute(vtkObject *object, unsigned long event, void *data) {
    this->Cursor->HandleEvent(object, event, data); };

protected:
  vtkSurfaceCursorCommand(vtkSurfaceCursor *cursor) {
    this->Cursor = cursor; };

  vtkSurfaceCursor* Cursor;

private:
  static vtkSurfaceCursorCommand *New(); // Not implemented.
  vtkSurfaceCursorCommand(); // Not implemented.
  vtkSurfaceCursorCommand(const vtkSurfaceCursorCommand&);  // Not implemented.
  void operator=(const vtkSurfaceCursorCommand&);  // Not implemented.
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
  this->Shapes = vtkDataSetCollection::New();
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

  this->MakeDefaultShapes();
  this->SetShape(0);

  this->Command = vtkSurfaceCursorCommand::New(this);
}

//----------------------------------------------------------------------------
vtkSurfaceCursor::~vtkSurfaceCursor()
{  
  this->SetRenderer(0);

  if (this->Command)
    {
    this->Command->Delete();
    }
  if (this->Matrix)
    {
    this->Matrix->Delete();
    }
  if (this->Shapes)
    {
    this->Shapes->Delete();
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
    this->Renderer->RemoveObserver(this->Command);
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
                                this->Command, -1);
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
  if ( (modifier & VTK_SCURSOR_SHIFT) )
    {
    return 1;
    }
  else if ( (modifier & VTK_SCURSOR_CONTROL) )
    {
    return 2;
    }

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
      if ((pickFlags & VTK_SCURSOR_IMAGE_ACTOR) ||
          ((pickFlags & VTK_SCURSOR_VOLUME) &&
           (pickFlags & (VTK_SCURSOR_CROP_PLANE | VTK_SCURSOR_CLIP_PLANE))))
        {
        shape = VTK_SCURSOR_CROSS_SPLIT;
        }
      else if ((pickFlags & VTK_SCURSOR_VOLUME) ||
               (pickFlags & VTK_SCURSOR_ACTOR))
        {
        shape = VTK_SCURSOR_CONE;
        }
      }
      break;

    case 1:
      {
      if ((pickFlags & (VTK_SCURSOR_CLIP_PLANE | VTK_SCURSOR_CROP_PLANE)))
        {
        shape = VTK_SCURSOR_PUSHER;
        }
      else if ((pickFlags & VTK_SCURSOR_PROP3D))
        {
        shape = VTK_SCURSOR_MOVER;
        }
      else
        {
        shape = VTK_SCURSOR_CROSSHAIRS;
        }
      }
      break;

    case 2:
      {
      if ((pickFlags & VTK_SCURSOR_CROP_PLANE))
        {
        shape = VTK_SCURSOR_ROCKER;
        }
      else if ((pickFlags & VTK_SCURSOR_PROP3D))
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
    case 1:
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

    case 2:
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
                                this->Mapper, this->Renderer);

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
void vtkSurfaceCursor::BindInteractor(vtkRenderWindowInteractor *iren)
{
  iren->AddObserver(vtkCommand::EnterEvent, this->Command);
  iren->AddObserver(vtkCommand::LeaveEvent, this->Command);
  iren->AddObserver(vtkCommand::KeyPressEvent, this->Command);
  iren->AddObserver(vtkCommand::KeyReleaseEvent, this->Command);
  iren->AddObserver(vtkCommand::MouseMoveEvent, this->Command);
  iren->AddObserver(vtkCommand::LeftButtonPressEvent, this->Command);
  iren->AddObserver(vtkCommand::RightButtonPressEvent, this->Command);
  iren->AddObserver(vtkCommand::MiddleButtonPressEvent, this->Command);
  iren->AddObserver(vtkCommand::LeftButtonReleaseEvent, this->Command);
  iren->AddObserver(vtkCommand::RightButtonReleaseEvent, this->Command);
  iren->AddObserver(vtkCommand::MiddleButtonReleaseEvent, this->Command);

  vtkRenderWindow *renwin = iren->GetRenderWindow();
  if (renwin)
    {
    renwin->AddObserver(vtkCommand::StartEvent, this->Command);
    }
  else
    {
    vtkErrorMacro("Connect the RenderWindow to the Interactor before calling"
                  " BindInteractor()");
    }
}

//----------------------------------------------------------------------------
void vtkSurfaceCursor::HandleEvent(vtkObject *object, unsigned long event,
                                   void *)
{
  // Capture the start of the RenderWindow and Renderer renders.
  if (event == vtkCommand::StartEvent)
    {
    if (object->IsA("vtkRenderWindow"))
      {
      // Set the DisplayPosition when the RenderWindow renders.  We will
      // only be hooked into RenderWindow::StartEvent if BindInteractor()
      // was called with a valid interactor.
      vtkRenderWindowInteractor *iren =
        static_cast<vtkRenderWindow *>(object)->GetInteractor();

      if (iren)
        {
        int x, y;
        iren->GetEventPosition(x, y);
        this->SetDisplayPosition(x, y);
        }
      }
    else if (object == this->Renderer)
      {
      // Compute the position when the Renderer renders, since it needs to
      // update all of the props in the scene.
      this->ComputePosition();
      }

    return;
    }

  // Check to see if the object is an interactor.
  vtkRenderWindowInteractor *iren =
    vtkRenderWindowInteractor::SafeDownCast(object);

  if (!iren)
    {
    return;
    } 

  // The interactor events are used to do just three things:
  // 1) set this->MouseInRenderer to control cursor visibility
  // 2) set this->Modifier for modifier keys and mouse buttons
  // 3) call MoveToDisplayPosition(x,y) when the mouse moves

  switch (event)
    {
    case vtkCommand::EnterEvent:
      {
      this->SetMouseInRenderer(1);
      iren->GetRenderWindow()->HideCursor();
      }
      break;
    case vtkCommand::LeaveEvent:
      {
      this->SetMouseInRenderer(0);
      iren->GetRenderWindow()->ShowCursor();
      }
      break;
    case vtkCommand::KeyPressEvent:
      {
      int modifierBits = this->ModifierFromKeySym(iren->GetKeySym());
      if (modifierBits)
        {
        this->SetModifier(this->Modifier | modifierBits);
        }
      }
      break;
    case vtkCommand::KeyReleaseEvent:
      {
      int modifierBits = this->ModifierFromKeySym(iren->GetKeySym());
      if (modifierBits)
        {
        this->SetModifier(this->Modifier & ~modifierBits);
        }
      }
      break;
    case vtkCommand::MouseMoveEvent:
      {
      int x, y;
      iren->GetEventPosition(x, y);
      // The Enter/Leave events aren't enough, because mouse drags don't
      // post the Leave event until the mouse button is released.
      int inRenderer = (this->Renderer && this->Renderer->IsInViewport(x, y));
      if (inRenderer != this->MouseInRenderer)
        {
        this->SetMouseInRenderer(inRenderer);
        if (inRenderer)
          {
          iren->GetRenderWindow()->HideCursor();
          }
        else
          {
          iren->GetRenderWindow()->ShowCursor();
          }
        }
      this->MoveToDisplayPosition(x, y);
      }
      break;
    case vtkCommand::LeftButtonPressEvent:
      {
      this->SetModifier(this->Modifier | VTK_SCURSOR_B1);
      }
      break;
    case vtkCommand::RightButtonPressEvent:
      {
      this->SetModifier(this->Modifier | VTK_SCURSOR_B2);
      }
      break;
    case vtkCommand::MiddleButtonPressEvent:
      {
      this->SetModifier(this->Modifier | VTK_SCURSOR_B3);
      }
      break;
    case vtkCommand::LeftButtonReleaseEvent:
      {
      this->SetModifier(this->Modifier & ~VTK_SCURSOR_B1);
      }
      break;
    case vtkCommand::RightButtonReleaseEvent:
      {
      this->SetModifier(this->Modifier & ~VTK_SCURSOR_B2);
      }
      break;
    case vtkCommand::MiddleButtonReleaseEvent:
      {
      this->SetModifier(this->Modifier & ~VTK_SCURSOR_B3);
      }
      break;
    }

  iren->Render();
}

//----------------------------------------------------------------------------
void vtkSurfaceCursor::MoveToDisplayPosition(double x, double y)
{
  switch (this->Action)
    {
    case VTK_SCURSOR_PUSH:
      {
      }
      break;
    case VTK_SCURSOR_ROTATE:
      {
      }
      break;
    case VTK_SCURSOR_SPIN:
      {
      }
      break;
    }

  this->SetDisplayPosition(x, y);
}

//----------------------------------------------------------------------------
void vtkSurfaceCursor::SetMouseInRenderer(int inside)
{
  if (this->MouseInRenderer == inside)
    {
    return;
    }
    
  this->MouseInRenderer = inside;
  this->Actor->SetVisibility(inside);
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
    // Terminate the previous action.
    }

  this->Action = action;
  this->Modified();

  if (action)
    {
    // Initiate the new action.
    }
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
  vtkDataSet *data = this->Shapes->GetItem(shape);

  if (data)
    {
    this->Mapper->SetInput(data);
    this->Shape = shape;
    }
}

//----------------------------------------------------------------------------
int vtkSurfaceCursor::AddShape(vtkDataSet *data)
{
  this->Shapes->AddItem(data);
  
  return (this->Shapes->GetNumberOfItems() - 1);
}

//----------------------------------------------------------------------------
void vtkSurfaceCursor::MakeDefaultShapes()
{
  vtkDataSet *data;

  data = this->MakePointerShape();
  this->AddShape(data);
  data->Delete();

  data = this->MakeCrosshairsShape();
  this->AddShape(data);
  data->Delete();

  data = this->MakeCrossShape(0);
  this->AddShape(data);
  data->Delete();

  data = this->MakeCrossShape(1);
  this->AddShape(data);
  data->Delete();

  data = this->MakeConeShape(0);
  this->AddShape(data);
  data->Delete();

  data = this->MakeConeShape(1);
  this->AddShape(data);
  data->Delete();

  data = this->MakeSphereShape(0);
  this->AddShape(data);
  data->Delete();

  data = this->MakeSphereShape(1);
  this->AddShape(data);
  data->Delete();

  data = this->MakeMoverShape(0);
  this->AddShape(data);
  data->Delete();

  data = this->MakeMoverShape(1);
  this->AddShape(data);
  data->Delete();

  data = this->MakePusherShape();
  this->AddShape(data);
  data->Delete();

  data = this->MakeSpinnerShape();
  this->AddShape(data);
  data->Delete();
}

//----------------------------------------------------------------------------
vtkDataSet *vtkSurfaceCursor::MakePointerShape()
{
  vtkUnsignedCharArray *scalars = vtkUnsignedCharArray::New();
  scalars->SetNumberOfComponents(4);
  vtkPoints *points = vtkPoints::New();
  vtkCellArray *strips = vtkCellArray::New();  
  vtkCellArray *lines = vtkCellArray::New();  
  
  static unsigned char black[4] = {  0,   0,   0, 255};
  static unsigned char white[4] = {255, 255, 255, 255};

  static double hotspot[2] = { 0.5, -0.5 };
  static double coords[7][2] = {
    {  1,  -1 },
    {  1, -15 },
    {  4, -12 },
    {  8, -19 },
    { 10, -18 },
    {  7, -11 },
    { 11, -11 },
  };

  static vtkIdType stripIds[] = {
    3, 0, 1, 2,
    3, 0, 2, 5,
    3, 0, 5, 6,
    4, 2, 3, 5, 4,
  };

  static vtkIdType lineIds[] = {
    8, 7, 8, 9, 10, 11, 12, 13, 7,
    8, 14, 15, 16, 17, 18, 19, 20, 14,
  };

  // Add the points three times: white, then black, and black again
  for (int i = 0; i < 7; i++)
    {
    points->InsertNextPoint(coords[i][0] - hotspot[0],
                            coords[i][1] - hotspot[1], 0.0);
    scalars->InsertNextTupleValue(white);
    }

  for (int j = 0; j < 7; j++)
    {
    points->InsertNextPoint(coords[j][0] - hotspot[0],
                            coords[j][1] - hotspot[1], +0.1);
    scalars->InsertNextTupleValue(black);
    }

  for (int k = 0; k < 7; k++)
    {
    points->InsertNextPoint(coords[k][0] - hotspot[0],
                            coords[k][1] - hotspot[1], -0.1);
    scalars->InsertNextTupleValue(black);
    }

  // Make the strips
  strips->InsertNextCell(stripIds[0], &stripIds[1]);
  strips->InsertNextCell(stripIds[4], &stripIds[5]);
  strips->InsertNextCell(stripIds[8], &stripIds[9]);
  strips->InsertNextCell(stripIds[12], &stripIds[13]);

  // Make the lines
  lines->InsertNextCell(lineIds[0], &lineIds[1]);
  lines->InsertNextCell(lineIds[8], &lineIds[9]);

  vtkPolyData *data = vtkPolyData::New();
  data->SetPoints(points);
  points->Delete();
  data->SetStrips(strips);
  strips->Delete();
  data->SetLines(lines);
  lines->Delete();
  data->GetPointData()->SetScalars(scalars);
  scalars->Delete();

  return data;
}

//----------------------------------------------------------------------------
vtkDataSet *vtkSurfaceCursor::MakeCrosshairsShape()
{
  vtkUnsignedCharArray *scalars = vtkUnsignedCharArray::New();
  scalars->SetNumberOfComponents(4);
  vtkPoints *points = vtkPoints::New();
  vtkCellArray *strips = vtkCellArray::New();  
  vtkCellArray *lines = vtkCellArray::New();  
  
  static unsigned char black[4] = {  0,   0,   0, 255};
  static unsigned char white[4] = {255, 255, 255, 255};

  const double radius = 8;
  const double inner = 3;

  static double coords[8][2] = {
    { 0, -radius }, { 0, -inner }, { 0, +inner }, { 0, +radius },
    { -radius, 0 }, { -inner, 0 }, { +inner, 0 }, { +radius, 0 },
  };

  static double outCoords[16][2] = {
    { -1, -radius-1 }, { +1, -radius-1 }, { +1, -inner+1 }, { -1, -inner+1 }, 
    { +1, +radius+1 }, { -1, +radius+1 }, { -1, +inner-1 }, { +1, +inner-1 }, 
    { -radius-1, +1 }, { -radius-1, -1 }, { -inner+1, -1 }, { -inner+1, +1 }, 
    { +radius+1, -1 }, { +radius+1, +1 }, { +inner-1, +1 }, { +inner-1, -1 }, 
  };

  static vtkIdType toplineIds[] = {
    2, 0, 1,
    2, 2, 3,
    2, 4, 5,
    2, 6, 7,
  };

  static vtkIdType botlineIds[] = {
    2, 8, 9,
    2, 10, 11,
    2, 12, 13,
    2, 14, 15,
  };

  static vtkIdType outlineIds[] = {
    5, 16, 17, 18, 19, 16,
    5, 20, 21, 22, 23, 20,
    5, 24, 25, 26, 27, 24,
    5, 28, 29, 30, 31, 28,
  };

  static vtkIdType stripIds[] = {
    4, 16, 17, 19, 18,
    4, 20, 21, 23, 22,
    4, 24, 25, 27, 26,
    4, 28, 29, 31, 30,
  };

  for (int i = 0; i < 8; i++)
    {
    points->InsertNextPoint(coords[i][0]+0.5, coords[i][1]-0.5, +0.1);
    scalars->InsertNextTupleValue(black);
    }

  for (int j = 0; j < 8; j++)
    {
    points->InsertNextPoint(coords[j][0]+0.5, coords[j][1]-0.5, -0.1);
    scalars->InsertNextTupleValue(black);
    }

  for (int k = 0; k < 16; k++)
    {
    points->InsertNextPoint(outCoords[k][0]+0.5, outCoords[k][1]-0.5, 0);
    scalars->InsertNextTupleValue(white);
    }

  // Make the crosshairs
  lines->InsertNextCell(toplineIds[0], &toplineIds[1]);
  lines->InsertNextCell(toplineIds[3], &toplineIds[4]);
  lines->InsertNextCell(toplineIds[6], &toplineIds[7]);
  lines->InsertNextCell(toplineIds[9], &toplineIds[10]);

  // Make the crosshairs
  lines->InsertNextCell(botlineIds[0], &botlineIds[1]);
  lines->InsertNextCell(botlineIds[3], &botlineIds[4]);
  lines->InsertNextCell(botlineIds[6], &botlineIds[7]);
  lines->InsertNextCell(botlineIds[9], &botlineIds[10]);

  // Make the outline
  lines->InsertNextCell(outlineIds[0], &outlineIds[1]);
  lines->InsertNextCell(outlineIds[6], &outlineIds[7]);
  lines->InsertNextCell(outlineIds[12], &outlineIds[13]);
  lines->InsertNextCell(outlineIds[18], &outlineIds[19]);
  lines->InsertNextCell(outlineIds[24], &outlineIds[25]);

  // Fill the outline
  strips->InsertNextCell(stripIds[0], &stripIds[1]);
  strips->InsertNextCell(stripIds[5], &stripIds[6]);
  strips->InsertNextCell(stripIds[10], &stripIds[11]);
  strips->InsertNextCell(stripIds[15], &stripIds[16]);

  vtkPolyData *data = vtkPolyData::New();
  data->SetPoints(points);
  points->Delete();
  data->SetLines(lines);
  lines->Delete();
  data->SetStrips(strips);
  strips->Delete();
  data->GetPointData()->SetScalars(scalars);
  scalars->Delete();

  return data;
}

//----------------------------------------------------------------------------
vtkDataSet *vtkSurfaceCursor::MakeSphereShape(int dual)
{
  double pi = vtkMath::DoublePi();
  double radius = 5.0;
  int resolution = 9;

  vtkIdType *pointIds = new vtkIdType[4*(resolution+1)];

  vtkIntArray *scalars = vtkIntArray::New();
  vtkDoubleArray *normals = vtkDoubleArray::New();
  normals->SetNumberOfComponents(3);
  vtkPoints *points = vtkPoints::New();
  vtkCellArray *strips = vtkCellArray::New();
  vtkIdType nPoints = 0;
 
  int colorIndex = 0;

  for (int i = 0; i < 2; i++)
    {
    // The sign (i.e. for top or bottom) is stored in s
    double s = 1 - 2*colorIndex;

    // The unit position vector of the point is stored in v
    double v[3];
    v[0] = 0;
    v[1] = 0;
    v[2] = s;

    points->InsertNextPoint(radius*v[0], radius*v[1], radius*v[2]);
    normals->InsertNextTupleValue(v);
    scalars->InsertNextTupleValue(&colorIndex);
      
    int n = (resolution + 1)/2;
    int m = 2*resolution;

    for (int j = 1; j <= n; j++)
      {
      double phi = pi*j/resolution;
      double r = sin(phi);
      v[2] = cos(phi)*s;
      if (2*j >= resolution)
        {
        v[2] = 0;
        }

      for (int k = 0; k < m; k++)
        {
        double theta = pi*k/resolution;
        v[0] = r*cos(theta);
        v[1] = r*sin(theta)*s;
        points->InsertNextPoint(radius*v[0], radius*v[1], radius*v[2]);
        normals->InsertNextTupleValue(v);
        scalars->InsertNextTupleValue(&colorIndex);
        }
      }

    // Make the fan for the top
    pointIds[0] = nPoints++;
    for (int ii = 0; ii < (m-1); ii++)
      {
      pointIds[1] = nPoints + ii;
      pointIds[2] = nPoints + ii + 1;
      strips->InsertNextCell(3, pointIds);
      }
    pointIds[1] = nPoints + m - 1;
    pointIds[2] = nPoints;
    strips->InsertNextCell(3, pointIds);

    // Make the strips for the rest
    for (int jj = 1; jj < n; jj++)
      {
      for (int kk = 0; kk < m; kk++)
        {
        pointIds[2*kk] = nPoints + kk;
        pointIds[2*kk+1] = nPoints + kk + m;
        }
      pointIds[2*m] = nPoints;
      pointIds[2*m+1] = nPoints + m;
      strips->InsertNextCell(2*(m+1), pointIds);
      nPoints += m;
      }

    nPoints += m;

    if (dual)
      {
      colorIndex = 1;
      }
    }

  delete [] pointIds;

  vtkPolyData *data = vtkPolyData::New();
  data->SetPoints(points);
  points->Delete();
  data->SetStrips(strips);
  strips->Delete();
  data->GetPointData()->SetScalars(scalars);
  scalars->Delete();
  data->GetPointData()->SetNormals(normals);
  normals->Delete();

  return data;
}

//----------------------------------------------------------------------------
vtkDataSet *vtkSurfaceCursor::MakeConeShape(int dual)
{
  double pi = vtkMath::DoublePi();
  double radius = 8.0;
  double height = 15.0;
  int resolution = 20;

  vtkIdType *pointIds = new vtkIdType[2*(resolution+1)];

  vtkIntArray *scalars = vtkIntArray::New();
  vtkDoubleArray *normals = vtkDoubleArray::New();
  normals->SetNumberOfComponents(3);
  vtkPoints *points = vtkPoints::New();
  vtkCellArray *strips = vtkCellArray::New();
  vtkIdType nPoints = 0;
 
  int sides = (dual ? 2 : 1);

  for (int colorIndex = 0; colorIndex < sides; colorIndex++)
    {
    // The sign (i.e. for top or bottom) is stored in s
    double s = 1 - 2*colorIndex;

    // The length of the side of the cone
    double l = sqrt(radius*radius + height*height);
    double f1 = radius/l;
    double f2 = height/l;

    // The unit normal vector
    double v[3];

    // The point of the cone
    for (int i = 0; i < 2; i++)
      {
      double r = radius*i;
      double z = height*i;
      double offset = 0.5*(1 - i);

      for (int j = 0; j < resolution; j++)
        {
        double theta = 2*pi*(j + offset)/resolution;
        double ct = cos(theta);
        double st = sin(theta);
        v[0] = f2*ct;
        v[1] = f2*st*s;
        v[2] = -f1*s;
        points->InsertNextPoint(r*ct, r*st*s, z*s);
        normals->InsertNextTupleValue(v);
        scalars->InsertNextTupleValue(&colorIndex);
        }
      }

    // The base of the cone
    v[0] = 0;
    v[1] = 0;
    v[2] = s;
    points->InsertNextPoint(0, 0, height*s);
    normals->InsertNextTupleValue(v);
    scalars->InsertNextTupleValue(&colorIndex); 

    for (int k = 0; k < resolution; k++)
      {
      double theta = 2*pi*k/resolution;
      points->InsertNextPoint(radius*cos(theta), radius*sin(theta)*s, height*s);
      normals->InsertNextTupleValue(v);
      scalars->InsertNextTupleValue(&colorIndex);
      }

    // Make the fan for the top
    for (int ii = 0; ii < (resolution-1); ii++)
      {
      pointIds[0] = nPoints + ii;
      pointIds[1] = nPoints + ii + resolution + 1;
      pointIds[2] = nPoints + ii + resolution;
      strips->InsertNextCell(3, pointIds);
      }
    pointIds[0] = nPoints + 2*resolution - 1;
    pointIds[1] = nPoints + resolution - 1;
    pointIds[2] = nPoints + resolution;
    strips->InsertNextCell(3, pointIds);
    nPoints += 2*resolution;
 
    // Make the fan for the base
    pointIds[0] = nPoints++;
    for (int jj = 0; jj < (resolution-1); jj++)
      {
      pointIds[1] = nPoints + jj;
      pointIds[2] = nPoints + jj + 1;
      strips->InsertNextCell(3, pointIds);
      }
    pointIds[1] = nPoints + resolution - 1;
    pointIds[2] = nPoints;
    strips->InsertNextCell(3, pointIds);
    nPoints += resolution;
    }

  delete [] pointIds;

  vtkPolyData *data = vtkPolyData::New();
  data->SetPoints(points);
  points->Delete();
  data->SetStrips(strips);
  strips->Delete();
  data->GetPointData()->SetScalars(scalars);
  scalars->Delete();
  data->GetPointData()->SetNormals(normals);
  normals->Delete();

  return data;
}

//----------------------------------------------------------------------------
vtkPolyData *vtkSurfaceCursor::MakeWarpedArrow(double warpX, double warpY,
                                               double warpZ, double warpScale)
{
  vtkPoints *points = vtkPoints::New();
  vtkCellArray *polys = vtkCellArray::New();

  static double coords[7][3] = {
    { 4, 2, 0 }, { 15, 2, 0 }, { 14, 8, 0 }, { 24, 0.01, 0 },
    { 14, -8, 0 }, { 15, -2, 0 }, { 4, -2, 0 },
  };

  static vtkIdType polyIds[] = {
    7, 0, 1, 2, 3, 4, 5, 6,
  };

  for (int i = 0; i < 7; i++)
    {
    points->InsertNextPoint(coords[i]);
    }
  polys->InsertNextCell(polyIds[0], &polyIds[1]);

  vtkPolyData *arrow = vtkPolyData::New();
  arrow->SetPoints(points);
  points->Delete();
  arrow->SetPolys(polys);
  polys->Delete();

  vtkImplicitModeller *modeller = vtkImplicitModeller::New();
  modeller->SetInput(arrow);
  modeller->SetSampleDimensions(32, 16, 8);

  vtkContourFilter *contour = vtkContourFilter::New();
  contour->SetInputConnection(modeller->GetOutputPort());
  contour->SetValue(0, 0.5);

  // The image is inside-out, and so is the contour
  vtkReverseSense *reverse = vtkReverseSense::New();
  reverse->SetInputConnection(contour->GetOutputPort());
  reverse->ReverseCellsOn();

  vtkWarpTo *warp = vtkWarpTo::New();
  warp->SetInputConnection(reverse->GetOutputPort());
  warp->SetPosition(warpX, warpY, warpZ);
  warp->SetScaleFactor(warpScale);
  warp->AbsoluteOn();

  vtkPolyDataNormals *polyNormals = vtkPolyDataNormals::New();
  polyNormals->SetInputConnection(warp->GetOutputPort());
  polyNormals->SplittingOn();

  vtkStripper *stripper = vtkStripper::New();
  stripper->SetInputConnection(polyNormals->GetOutputPort());

  stripper->Update();

  vtkPolyData *data = vtkPolyData::New();
  data->DeepCopy(stripper->GetOutput());

  polyNormals->Delete();
  stripper->Delete();
  warp->Delete();
  reverse->Delete();
  contour->Delete();
  modeller->Delete();

  return data;
}

//----------------------------------------------------------------------------
vtkDataSet *vtkSurfaceCursor::MakeMoverShape(int warped)
{
  vtkPolyData *leafData = vtkSurfaceCursor::MakeWarpedArrow(
    10, 0, -10, (warped ? 1.0 : 0.5));
  vtkPoints *leafPoints = leafData->GetPoints();
  vtkCellArray *leafStrips = leafData->GetStrips();
  vtkDataArray *leafNormals = leafData->GetPointData()->GetNormals();

  vtkPolyData *data = vtkPolyData::New();
  vtkPoints *points = vtkPoints::New();
  vtkDoubleArray *normals = vtkDoubleArray::New();
  normals->SetNumberOfComponents(3);
  vtkCellArray *strips = vtkCellArray::New();
  vtkIntArray *scalars = vtkIntArray::New();

  vtkTransform *transform = vtkTransform::New();
  transform->PostMultiply();

  static double rotate90[16] = {
    0,-1, 0, 0,
    1, 0, 0, 0,
    0, 0, 1, 0,
    0, 0, 0, 1,
  };

  if (warped)
    {
    transform->RotateY(36);
    transform->Scale(1.5, 1.0, 1.0);
    transform->Translate(0, 0, 14);
    }

  int color = 0;

  for (int j = 0; j < 4; j++)
    {
    vtkIdType nPoints = points->GetNumberOfPoints();
    transform->TransformPoints(leafPoints, points);
    transform->TransformNormals(leafNormals, normals);
    vtkIdType npts;
    vtkIdType *pts;
    leafStrips->InitTraversal();
    while (leafStrips->GetNextCell(npts, pts))
      {
      strips->InsertNextCell(npts);
      for (vtkIdType k = 0; k < npts; k++)
        {
        strips->InsertCellPoint(pts[k] + nPoints);
        }
      }
    vtkIdType n = leafPoints->GetNumberOfPoints();
    for (int ii = 0; ii < n; ii++)
      {
      scalars->InsertNextTupleValue(&color);
      }
    transform->Concatenate(rotate90);
    color = !color;
    }
  transform->Delete();

  leafData->Delete();

  data->SetPoints(points);
  points->Delete();
  data->GetPointData()->SetNormals(normals);
  normals->Delete();
  data->GetPointData()->SetScalars(scalars);
  scalars->Delete();
  data->SetStrips(strips);
  strips->Delete();

  return data;
}

//----------------------------------------------------------------------------
vtkDataSet *vtkSurfaceCursor::MakePusherShape()
{
  vtkPolyData *leafData = vtkSurfaceCursor::MakeWarpedArrow(
    10.0, 0.0, 10.0, 0.0);
  vtkPoints *leafPoints = leafData->GetPoints();
  vtkCellArray *leafStrips = leafData->GetStrips();
  vtkDataArray *leafNormals = leafData->GetPointData()->GetNormals();

  vtkPolyData *data = vtkPolyData::New();
  vtkPoints *points = vtkPoints::New();
  vtkDoubleArray *normals = vtkDoubleArray::New();
  normals->SetNumberOfComponents(3);
  vtkCellArray *strips = vtkCellArray::New();
  vtkCellArray *lines = vtkCellArray::New();
  vtkIntArray *scalars = vtkIntArray::New();

  vtkTransform *transform = vtkTransform::New();
  transform->PostMultiply();

  static double rotate90[16] = {
    0, 0, 1, 0,
    0, 1, 0, 0,
   -1, 0, 0, 0,
    0, 0, 0, 1,
  };

  transform->Concatenate(rotate90);
  transform->Scale(1.0, 1.0, 1.0);
  transform->Translate(0, 0, 24);

  int color = 0;

  for (int i = 0; i < 2; i++)
    {
    vtkIdType nPoints = points->GetNumberOfPoints();
    transform->TransformPoints(leafPoints, points);
    transform->TransformNormals(leafNormals, normals);
    vtkIdType npts;
    vtkIdType *pts;
    leafStrips->InitTraversal();
    while (leafStrips->GetNextCell(npts, pts))
      {
      strips->InsertNextCell(npts);
      for (vtkIdType j = 0; j < npts; j++)
        {
        strips->InsertCellPoint(pts[j] + nPoints);
        }
      }
    vtkIdType nn = leafPoints->GetNumberOfPoints();
    for (int ii = 0; ii < nn; ii++)
      {
      scalars->InsertNextTupleValue(&color);
      }
    // This transform puts the arrows tail-to-tail, instead of tip-to-tip 
    transform->Translate(0, 0, -44);
    transform->Scale(-1,1,-1);
    color = !color;
    }
  transform->Delete();

  leafData->Delete();

  // Make a ring for when arrow viewed tail-on
  const double lineRadius = 8;
  const int lineResolution = 24; // must be divisible by 8
  int polylen = lineResolution/8 + 1;
  double normal[3];
  normal[0] = normal[1] = 0.0;
  normal[2] = 1.0;
  for (int j = 0; j < 2; j++)
    {
    color = 0;
    for (int k = 0; k < 8; k++)
      {
      vtkIdType nPoints = points->GetNumberOfPoints();
      lines->InsertNextCell(polylen);
      for (int ii = 0; ii < polylen; ii++)
        {
        double angle = 2*vtkMath::DoublePi()*(k*(polylen-1)+ii)/lineResolution;
        // The 0.99 is to make x slightly thinner than y
        points->InsertNextPoint(0.99*lineRadius*cos(angle),
                                lineRadius*sin(angle),
                                0.1*(1 - 2*j));
        scalars->InsertNextTupleValue(&color);
        normals->InsertNextTupleValue(normal);
        lines->InsertCellPoint(ii + nPoints);
        }
      normal[2] = -normal[2];
      color = !color;
      }
    }

  data->SetPoints(points);
  points->Delete();
  data->GetPointData()->SetNormals(normals);
  normals->Delete();
  data->GetPointData()->SetScalars(scalars);
  scalars->Delete();
  data->SetStrips(strips);
  strips->Delete();
  data->SetLines(lines);
  lines->Delete();

  return data;
}

//----------------------------------------------------------------------------
vtkDataSet *vtkSurfaceCursor::MakeSpinnerShape()
{
  vtkPolyData *leafData = vtkSurfaceCursor::MakeWarpedArrow(
    10.0, 0.0, -5.0, 1.0);
  vtkPoints *leafPoints = leafData->GetPoints();
  vtkCellArray *leafStrips = leafData->GetStrips();
  vtkDataArray *leafNormals = leafData->GetPointData()->GetNormals();

  vtkPolyData *data = vtkPolyData::New();
  vtkPoints *points = vtkPoints::New();
  vtkDoubleArray *normals = vtkDoubleArray::New();
  normals->SetNumberOfComponents(3);
  vtkCellArray *strips = vtkCellArray::New();
  vtkIntArray *scalars = vtkIntArray::New();

  vtkTransform *transform = vtkTransform::New();
  transform->PostMultiply();

  static double rotate90[16] = {
    1, 0, 0, 0,
    0, 0, 1, 0,
    0,-1, 0, 0,
    0, 0, 0, 1,
  };

  transform->Concatenate(rotate90);
  transform->Scale(3.0, 3.0, 1.5);
  transform->Translate(-30, 16, 3.75);

  int color = 0;

  for (int j = 0; j < 2; j++)
    {
    vtkIdType nPoints = points->GetNumberOfPoints();
    transform->TransformPoints(leafPoints, points);
    transform->TransformNormals(leafNormals, normals);
    vtkIdType npts;
    vtkIdType *pts;
    leafStrips->InitTraversal();
    while (leafStrips->GetNextCell(npts, pts))
      {
      strips->InsertNextCell(npts);
      for (vtkIdType k = 0; k < npts; k++)
        {
        strips->InsertCellPoint(pts[k] + nPoints);
        }
      }
    vtkIdType n = leafPoints->GetNumberOfPoints();
    for (int ii = 0; ii < n; ii++)
      {
      scalars->InsertNextTupleValue(&color);
      }
    transform->Scale(-1, -1, 1);
    color = !color;
    }
  transform->Delete();

  leafData->Delete();

  data->SetPoints(points);
  points->Delete();
  data->GetPointData()->SetNormals(normals);
  normals->Delete();
  data->GetPointData()->SetScalars(scalars);
  scalars->Delete();
  data->SetStrips(strips);
  strips->Delete();

  return data;
}

//----------------------------------------------------------------------------
vtkDataSet *vtkSurfaceCursor::MakeCrossShape(int dual)
{
  double radius = 10.0;
  double inner = 3.5;
  double thickness = 2.0;

  double xmin = inner;
  double xmax = radius;
  double ymin = -thickness*0.5;
  double ymax = +thickness*0.5;
  double zmin = 0;
  double zmax = thickness*0.5;

  vtkIntArray *scalars = vtkIntArray::New();
  vtkPoints *points = vtkPoints::New();
  vtkCellArray *strips = vtkCellArray::New();
  vtkIdType nPoints = 0;
  
  int colorIndex = 0;

  for (int i = 0; i < 2; i++)
    {
    for (int j = 0; j < 4; j++)
      {
      double z = zmax;
      for (int k = 0; k < 2; k++)
        {
        points->InsertNextPoint(xmin, ymin, z);
        points->InsertNextPoint(xmin, ymax, z);
        points->InsertNextPoint(xmax, ymax, z);
        points->InsertNextPoint(xmax, ymin, z);
        scalars->InsertNextTupleValue(&colorIndex);
        scalars->InsertNextTupleValue(&colorIndex);
        scalars->InsertNextTupleValue(&colorIndex);
        scalars->InsertNextTupleValue(&colorIndex);
        z = zmin;
        }
   
      static vtkIdType rectIds[5][4] = { { 1, 0, 2, 3 },
                                         { 4, 0, 5, 1 },
                                         { 5, 1, 6, 2 },
                                         { 6, 2, 7, 3 },
                                         { 7, 3, 4, 0 } };
      vtkIdType pointIds[4];
      for (int ii = 0; ii < 5; ii++)
        {
        for (int jj = 0; jj < 4; jj++)
          {
          pointIds[jj] = rectIds[ii][jj]+nPoints;
          }
        strips->InsertNextCell(4, pointIds);
        }
      nPoints += 8;

      // do a rotation of 90 degrees for next piece
      double tmp1 = ymin;
      double tmp2 = ymax;
      ymin = -xmax;
      ymax = -xmin;
      xmin = tmp1;
      xmax = tmp2;
      }

    // do the other side
    zmin = -zmin;
    zmax = -zmax;
    xmin = -xmin;
    xmax = -xmax;

    if (dual)
      {
      colorIndex = 1;
      }
    }

  vtkPolyData *data = vtkPolyData::New();
  data->SetPoints(points);
  points->Delete();
  data->SetStrips(strips);
  strips->Delete();
  data->GetPointData()->SetScalars(scalars);
  scalars->Delete();

  return data;
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
                                               vtkRenderer *renderer)
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
  double bounds[6], thickness[3];
  cursorMapper->GetBounds(bounds);
  thickness[2] = bounds[5] - bounds[4];
  int minDim = 2;
  for (int i = 0; i < 2; i++)
    {
    thickness[i] = bounds[2*i+1] - bounds[2*i];
    if (thickness[i] < thickness[minDim])
      {
      minDim = i;
      }
    }
  if (minDim != 2 && thickness[minDim] < 0.5*thickness[2])
    {
    direction = minDim;
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

//----------------------------------------------------------------------------
int vtkSurfaceCursor::ModifierFromKeySym(const char *keysym)
{
  if (keysym)
    {
    // These match the Tk modifier bits.  Also the following:
    // 1st button = 256, 2nd button = 512, middle button = 1024
    if (strncmp(keysym, "Shift_", 6) == 0)
      {
      return VTK_SCURSOR_SHIFT;
      }
    else if (strncmp(keysym, "Caps_Lock", 9) == 0)
      {
      return VTK_SCURSOR_CAPS;
      }
    else if (strncmp(keysym, "Control_", 8) == 0)
      {
      return VTK_SCURSOR_CONTROL;
      }
    else if (strncmp(keysym, "Meta_", 5) == 0)
      {
      return VTK_SCURSOR_META;
      }
    else if (strncmp(keysym, "Alt_", 4) == 0)
      {
      return VTK_SCURSOR_ALT;
      }
    }

  return 0;
}


