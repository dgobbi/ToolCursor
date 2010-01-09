/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkSurfaceCursorInteractorObserver.cxx,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkSurfaceCursorInteractorObserver.h"
#include "vtkObjectFactory.h"

#include "vtkSurfaceCursor.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkInteractorStyle.h"
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkCallbackCommand.h"

vtkCxxRevisionMacro(vtkSurfaceCursorInteractorObserver, "$Revision: 1.2 $");
vtkStandardNewMacro(vtkSurfaceCursorInteractorObserver);

vtkCxxSetObjectMacro(vtkSurfaceCursorInteractorObserver,SurfaceCursor, vtkSurfaceCursor);

//----------------------------------------------------------------------------
vtkSurfaceCursorInteractorObserver::vtkSurfaceCursorInteractorObserver()
{
  // The surface cursor that this object handles the events for
  this->SurfaceCursor = 0;

  // Set priority to be higher than the InteractorStyle
  this->Priority = 0.1;

  // InteractorObservers use a static method to handle window events.
  this->EventCallbackCommand->SetCallback(
    vtkSurfaceCursorInteractorObserver::ProcessEvents);

  // Make a callback for the RenderWindow render start & end
  this->PassiveEventCallbackCommand = vtkCallbackCommand::New();
  this->PassiveEventCallbackCommand->SetPassiveObserver(1);
  this->PassiveEventCallbackCommand->SetClientData(this);
  this->PassiveEventCallbackCommand->SetCallback(
    vtkSurfaceCursorInteractorObserver::ProcessPassiveEvents);
}

//----------------------------------------------------------------------------
vtkSurfaceCursorInteractorObserver::~vtkSurfaceCursorInteractorObserver()
{
  this->PassiveEventCallbackCommand->Delete();
  this->SurfaceCursor->Delete();
}

//----------------------------------------------------------------------------
void vtkSurfaceCursorInteractorObserver::PrintSelf(ostream& os,
                                                   vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "SurfaceCursor: " << this->SurfaceCursor << "\n";
}

//----------------------------------------------------------------------------
void vtkSurfaceCursorInteractorObserver::SetEnabled(int enable)
{
  vtkRenderWindowInteractor *iren = this->Interactor;
  if (!iren)
    {
    vtkErrorMacro("The interactor must be set prior to enabling/"
                  "disabling the widget");
    return;
    }

  vtkRenderWindow *renwin = iren->GetRenderWindow();
  if (!renwin)
    {
    vtkErrorMacro("Connect the RenderWindow to the Interactor before"
                  " enabling the widget");
    return;
    }

  if (enable && !this->Enabled)
    {
    this->Enabled = 1;

    vtkCommand *command = this->EventCallbackCommand;
    float priority = this->Priority;

    iren->AddObserver(vtkCommand::KeyPressEvent, command, priority);
    iren->AddObserver(vtkCommand::KeyReleaseEvent, command, priority);
    iren->AddObserver(vtkCommand::LeaveEvent, command, priority);
    iren->AddObserver(vtkCommand::MouseMoveEvent, command, priority);
    iren->AddObserver(vtkCommand::LeftButtonPressEvent, command, priority);
    iren->AddObserver(vtkCommand::RightButtonPressEvent, command, priority);
    iren->AddObserver(vtkCommand::MiddleButtonPressEvent, command, priority);
    iren->AddObserver(vtkCommand::LeftButtonReleaseEvent, command, priority);
    iren->AddObserver(vtkCommand::RightButtonReleaseEvent, command, priority);
    iren->AddObserver(vtkCommand::MiddleButtonReleaseEvent, command, priority);

    command = this->PassiveEventCallbackCommand;
    renwin->AddObserver(vtkCommand::StartEvent, command);
    renwin->AddObserver(vtkCommand::EndEvent, command);
    iren->AddObserver(vtkCommand::LeaveEvent, command);
    iren->AddObserver(vtkCommand::MouseMoveEvent, command);
    iren->AddObserver(vtkCommand::KeyPressEvent, command);
    iren->AddObserver(vtkCommand::KeyReleaseEvent, command);

    this->InvokeEvent(vtkCommand::EnableEvent, NULL);
    }
  else if (!enable && this->Enabled)
    {
    this->Enabled = 0;

    iren->RemoveObserver(this->EventCallbackCommand);
    iren->RemoveObserver(this->PassiveEventCallbackCommand);
    renwin->RemoveObserver(this->PassiveEventCallbackCommand);

    this->InvokeEvent(vtkCommand::DisableEvent, NULL);
    }
}

//----------------------------------------------------------------------------
void vtkSurfaceCursorInteractorObserver::ProcessPassiveEvents(
  vtkObject *object, unsigned long event, void *clientdata, void *)
{
  vtkSurfaceCursorInteractorObserver* self =
    reinterpret_cast<vtkSurfaceCursorInteractorObserver *>(clientdata);

  vtkSurfaceCursor *cursor = self->GetSurfaceCursor();
  vtkRenderWindowInteractor *iren = self->GetInteractor();

  // Look for events from the RenderWindow
  vtkRenderWindow *renwin = vtkRenderWindow::SafeDownCast(object);
  if (renwin && iren && renwin->GetInteractor() == iren)
    {
    if (event == vtkCommand::StartEvent)
      {
      // Just before the RenderWindow renders, get the mouse position and
      // use it to set the cursor position.

      int x, y;
      iren->GetEventPosition(x, y);
      cursor->SetDisplayPosition(x, y);
      }
    else if (event == vtkCommand::EndEvent)
      {
      // At end of RenderWindow render, check whether cursor is visible.
      // Hide system cursor if 3D cursor is visible.
      if (cursor->GetVisibility())
        {
        renwin->HideCursor();
        }
      else
        {
        renwin->ShowCursor();
        }
      }
   }

  // Look for events from the interactor
  else if (iren && vtkRenderWindowInteractor::SafeDownCast(object) == iren)
    {
    if (event == vtkCommand::MouseMoveEvent)
      {
      // Need to check if mouse is in the renderer, even if some other
      // observer has focus, so do it here as a passive operation.

      vtkRenderer *renderer = cursor->GetRenderer();
      int x, y;
      iren->GetEventPosition(x, y);
      int inRenderer = (renderer && renderer->IsInViewport(x, y));
      if (inRenderer != cursor->GetMouseInRenderer())
        {
        cursor->SetMouseInRenderer(inRenderer);
        }
      }
    else if (event == vtkCommand::LeaveEvent)
      {
      // Mouse move events might cease after mouse leaves the window,
      // leaving EventPosition with the last in-window value

      if (cursor->GetMouseInRenderer())
        {
        cursor->SetMouseInRenderer(0);
        }
      }
    else if (event == vtkCommand::KeyPressEvent)
      {
      // We need to know the exact moment when modifier keys change
      int modifierBits = self->ModifierFromKeySym(iren->GetKeySym());
      if (modifierBits)
        {
        cursor->SetModifier(cursor->GetModifier() | modifierBits);
        }
      }
    else if (event == vtkCommand::KeyReleaseEvent)
      {
      int modifierBits = self->ModifierFromKeySym(iren->GetKeySym());
      if (modifierBits)
        {
        cursor->SetModifier(cursor->GetModifier() & ~modifierBits);
        }
      }
    }
}

//----------------------------------------------------------------------------
void vtkSurfaceCursorInteractorObserver::ProcessEvents(vtkObject *object,
                                                       unsigned long event,
                                                       void *clientdata,
                                                       void *)
{
  vtkSurfaceCursorInteractorObserver* self =
    reinterpret_cast<vtkSurfaceCursorInteractorObserver *>(clientdata);

  vtkSurfaceCursor *cursor = self->GetSurfaceCursor();
  vtkRenderWindowInteractor *iren = self->GetInteractor();
  if (!iren || iren != vtkRenderWindowInteractor::SafeDownCast(object))
    {
    return;
    } 

  // Is it safe to grab the focus for the cursor?  Check to see if the
  // InteractorStyle is currently doing an action.
  vtkInteractorStyle *istyle =
    vtkInteractorStyle::SafeDownCast(iren->GetInteractorStyle());
  int allowTakeFocus = (istyle && istyle->GetState() == VTKIS_NONE);

  // The interactor events are used to do just three things:
  // 1) call cursor->SetMouseInRenderer() to control cursor visibility,
  //    see ProcessPassiveEvents
  // 2) call cursor->SetModifier() for modifier keys and mouse buttons,
  //    see ProcessPassiveEvents for keys
  // 3) call cursor->MoveToDisplayPosition(x,y) when the mouse moves

  switch (event)
    {
    case vtkCommand::MouseMoveEvent:
      {
      int x, y;
      iren->GetEventPosition(x, y);
      cursor->MoveToDisplayPosition(x, y);
      }
      break;
    case vtkCommand::LeftButtonPressEvent:
    case vtkCommand::RightButtonPressEvent:
    case vtkCommand::MiddleButtonPressEvent:
      {
      int button = 0;
      if (event == vtkCommand::LeftButtonPressEvent)
        {
        button = VTK_SCURSOR_B1;
        }
      else if (event == vtkCommand::RightButtonPressEvent)
        {
        button = VTK_SCURSOR_B2;
        }
      else if (event == vtkCommand::MiddleButtonPressEvent)
        {
        button = VTK_SCURSOR_B3;
        }

      if (allowTakeFocus && cursor->GetRequestingFocus())
        {
        self->GrabFocus(self->EventCallbackCommand,
                        self->EventCallbackCommand);
        // Make sure that no other observers see the event
        self->EventCallbackCommand->SetAbortFlag(1);

        int x, y;
        iren->GetEventPosition(x, y);
        cursor->SetDisplayPosition(x, y);
        cursor->SetModifier(cursor->GetModifier() & ~button);
        cursor->SetModifier(cursor->GetModifier() | button);
        }
      }
      break;
    case vtkCommand::LeftButtonReleaseEvent:
    case vtkCommand::RightButtonReleaseEvent:
    case vtkCommand::MiddleButtonReleaseEvent:
      {
      int button = 0;
      if (event == vtkCommand::LeftButtonReleaseEvent)
        {
        button = VTK_SCURSOR_B1;
        }
      else if (event == vtkCommand::RightButtonReleaseEvent)
        {
        button = VTK_SCURSOR_B2;
        }
      else if (event == vtkCommand::MiddleButtonReleaseEvent)
        {
        button = VTK_SCURSOR_B3;
        }

      if (allowTakeFocus)
        {
        self->ReleaseFocus();
        }

      int x, y;
      iren->GetEventPosition(x, y);
      cursor->SetDisplayPosition(x, y);
      cursor->SetModifier(cursor->GetModifier() & ~button);
      }
      break;
    case vtkCommand::LeaveEvent:
    case vtkCommand::KeyPressEvent:
    case vtkCommand::KeyReleaseEvent:
      // Just render for these events, they were handled as passive events
      break;
    default:
      return;
    }

  iren->Render();
}

//----------------------------------------------------------------------------
int vtkSurfaceCursorInteractorObserver::ModifierFromKeySym(const char *keysym)
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
