/*=========================================================================

  Program:   ToolCursor
  Module:    vtkWindowLevelTool.cxx

  Copyright (c) 2010 David Gobbi
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkWindowLevelTool.h"
#include "vtkObjectFactory.h"

#include "vtkToolCursor.h"
#include "vtkVolumePicker.h"
#include "vtkTransform.h"
#include "vtkRenderer.h"
#include "vtkImageProperty.h"
#include "vtkImageSlice.h"
#include "vtkCommand.h"

vtkStandardNewMacro(vtkWindowLevelTool);

//----------------------------------------------------------------------------
vtkWindowLevelTool::vtkWindowLevelTool()
{
  this->StartWindowLevel[0] = 1.0;
  this->StartWindowLevel[1] = 0.5;
}

//----------------------------------------------------------------------------
vtkWindowLevelTool::~vtkWindowLevelTool()
{
}

//----------------------------------------------------------------------------
void vtkWindowLevelTool::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

//----------------------------------------------------------------------------
void vtkWindowLevelTool::StartAction()
{
  this->Superclass::StartAction();

  if (this->CurrentImageProperty)
    {
    vtkImageProperty *property = this->CurrentImageProperty;
    this->StartWindowLevel[0] = property->GetColorWindow();
    this->StartWindowLevel[1] = property->GetColorLevel();
    }
}

//----------------------------------------------------------------------------
void vtkWindowLevelTool::StopAction()
{
  this->Superclass::StopAction();
  this->InvokeEvent(vtkCommand::EndWindowLevelEvent, this);
}

//----------------------------------------------------------------------------
void vtkWindowLevelTool::DoAction()
{
  this->Superclass::DoAction();

  vtkToolCursor *cursor = this->GetToolCursor();

  // Get the display position.
  double x, y, x0, y0;
  this->GetStartDisplayPosition(x0, y0);
  this->GetDisplayPosition(x, y);

  int *size = cursor->GetRenderer()->GetSize();

  double window = this->StartWindowLevel[0];
  double level = this->StartWindowLevel[1];

  level += window*(y - y0)*2.0/(size[1] + 1);
  window *= pow(10.0, (x - x0)*2.0/(size[1] + 1));

  if (this->CurrentImageProperty)
    {
    this->CurrentImageProperty->SetColorWindow(window);
    this->CurrentImageProperty->SetColorLevel(level);
    }
}
