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

vtkCxxRevisionMacro(vtkSurfaceCursorInteractorObserver, "$Revision: 1.1 $");
vtkStandardNewMacro(vtkSurfaceCursorInteractorObserver);

//----------------------------------------------------------------------------
vtkSurfaceCursorInteractorObserver::vtkSurfaceCursorInteractorObserver()
{
  this->SurfaceCursor = 0;

  // InteractorObservers use a static method to handle window events.
  // This callback will not receive keyboard events, the InteractorObserver
  // base class uses a separate callback for the keyboard.
  this->EventCallbackCommand->SetCallback(
    vtkSurfaceCursorInteractorObserver::ProcessEvents);
}

//----------------------------------------------------------------------------
vtkSurfaceCursorInteractorObserver::~vtkSurfaceCursorInteractorObserver()
{
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

  vtkCommand *command = this->EventCallbackCommand;
  float priority = this->Priority;

  if (enable && !this->Enabled)
    {
    this->Enabled = 1;

    iren->AddObserver(vtkCommand::EnterEvent, command, priority);
    iren->AddObserver(vtkCommand::LeaveEvent, command, priority);
    iren->AddObserver(vtkCommand::KeyPressEvent, command, priority);
    iren->AddObserver(vtkCommand::KeyReleaseEvent, command, priority);
    iren->AddObserver(vtkCommand::MouseMoveEvent, command, priority);
    iren->AddObserver(vtkCommand::LeftButtonPressEvent, command, priority);
    iren->AddObserver(vtkCommand::RightButtonPressEvent, command, priority);
    iren->AddObserver(vtkCommand::MiddleButtonPressEvent, command, priority);
    iren->AddObserver(vtkCommand::LeftButtonReleaseEvent, command, priority);
    iren->AddObserver(vtkCommand::RightButtonReleaseEvent, command, priority);
    iren->AddObserver(vtkCommand::MiddleButtonReleaseEvent, command, priority);

    renwin->AddObserver(vtkCommand::StartEvent, command);
    renwin->AddObserver(vtkCommand::EndEvent, command);

    this->InvokeEvent(vtkCommand::EnableEvent, NULL);
    }
  else if (!enable && this->Enabled)
    {
    this->Enabled = 0;

    iren->RemoveObserver(command);
    renwin->RemoveObserver(command);

    this->InvokeEvent(vtkCommand::DisableEvent, NULL);
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

  // Is it safe to grab the focus for the cursor?  Check to see if the
  // InteractorStyle is currently doing an action.
  int allowTakeFocus = 1;
  if (iren)
    {
    vtkInteractorStyle *istyle =
      vtkInteractorStyle::SafeDownCast(iren->GetInteractorStyle());
    allowTakeFocus = (istyle && istyle->GetState() == VTKIS_NONE);
    }

  // Look for events from the RenderWindow
  if (event == vtkCommand::StartEvent)
    {
    // Just before the RenderWindow renders, get the mouse position and
    // use it to set the cursor position.

    vtkRenderWindow *renwin = vtkRenderWindow::SafeDownCast(object);
    if (renwin)
      {
      vtkRenderWindowInteractor *iren = renwin->GetInteractor();
      if (iren && iren == self->GetInteractor())
        {
        int x, y;
        iren->GetEventPosition(x, y);
        cursor->SetDisplayPosition(x, y);
        }
      }

    // Must return to avoid the Render at the end of HandleEvent.
    return;
    }
  else if (event == vtkCommand::EndEvent)
    {
    // At end of RenderWindow render, check whether cursor is visible

    vtkRenderWindow *renwin = vtkRenderWindow::SafeDownCast(object);
    if (renwin)
      {
      // Hide system cursor if 3D cursor is visible.
      if (cursor->GetVisibility())
        {
        renwin->HideCursor();
        }
      else
        {
        renwin->ShowCursor();
        }

      // If something interesting is under the cursor, then take focus
      if (allowTakeFocus)
        {
        if (cursor->GetRequestingFocus())
          {
          self->GrabFocus(self->EventCallbackCommand,
                          self->EventCallbackCommand);
          }
        else
          {
          self->ReleaseFocus();
          }
        }
      }
    // Must return to avoid the Render at the end of HandleEvent.
    return;
    }

  // If the object is an interactor, then handle interaction events.
  if (!iren || iren != vtkRenderWindowInteractor::SafeDownCast(object))
    {
    return;
    } 

  // The interactor events are used to do just three things:
  // 1) call cursor->SetMouseInRenderer() to control cursor visibility
  // 2) call cursor->SetModifier() for modifier keys and mouse buttons
  // 3) call cursor->MoveToDisplayPosition(x,y) when the mouse moves

  switch (event)
    {
    case vtkCommand::EnterEvent:
      {
      cursor->SetMouseInRenderer(1);
      }
      break;
    case vtkCommand::LeaveEvent:
      {
      cursor->SetMouseInRenderer(0);
      }
      break;
    case vtkCommand::KeyPressEvent:
      {
      // Monitor the modifier keys like Shift, Ctrl, or Alt 
      int modifierBits = self->ModifierFromKeySym(iren->GetKeySym());
      if (modifierBits)
        {
        cursor->SetModifier(cursor->GetModifier() | modifierBits);
        }
      }
      break;
    case vtkCommand::KeyReleaseEvent:
      {
      int modifierBits = self->ModifierFromKeySym(iren->GetKeySym());
      if (modifierBits)
        {
        cursor->SetModifier(cursor->GetModifier() & ~modifierBits);
        }
      }
      break;
    case vtkCommand::MouseMoveEvent:
      {
      int x, y;
      iren->GetEventPosition(x, y);
      // The Enter/Leave events aren't enough, because mouse drags don't
      // post the Leave event until the mouse button is released.
      int inRenderer = (cursor->GetRenderer() &&
                        cursor->GetRenderer()->IsInViewport(x, y));

      if (inRenderer != cursor->GetMouseInRenderer())
        {
        cursor->SetMouseInRenderer(inRenderer);
        }
      // Note: the above checks should go in the cursor class
      cursor->MoveToDisplayPosition(x, y);
      }
      break;
    case vtkCommand::LeftButtonPressEvent:
      {
      cursor->SetModifier(cursor->GetModifier() | VTK_SCURSOR_B1);
      }
      break;
    case vtkCommand::RightButtonPressEvent:
      {
      cursor->SetModifier(cursor->GetModifier() | VTK_SCURSOR_B2);
      }
      break;
    case vtkCommand::MiddleButtonPressEvent:
      {
      cursor->SetModifier(cursor->GetModifier() | VTK_SCURSOR_B3);
      }
      break;
    case vtkCommand::LeftButtonReleaseEvent:
      {
      cursor->SetModifier(cursor->GetModifier() & ~VTK_SCURSOR_B1);
      }
      break;
    case vtkCommand::RightButtonReleaseEvent:
      {
      cursor->SetModifier(cursor->GetModifier() & ~VTK_SCURSOR_B2);
      }
      break;
    case vtkCommand::MiddleButtonReleaseEvent:
      {
      cursor->SetModifier(cursor->GetModifier() & ~VTK_SCURSOR_B3);
      }
      break;
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
