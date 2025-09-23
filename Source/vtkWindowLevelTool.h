/*=========================================================================

  Program:   ToolCursor
  Module:    vtkWindowLevelTool.h

  Copyright (c) 2010 David Gobbi
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkWindowLevelTool - Adjusts the Window/Level of an image.
// .SECTION Description
// This class adjusts the Window/Level of an image property.  The Window
// increases as the mouse moves right, the Level increases as the mouse
// moves up.

#ifndef vtkWindowLevelTool_h
#define vtkWindowLevelTool_h

#include "vtkToolCursorModule.h" // For export macro
#include "vtkTool.h"

class vtkImageProperty;

class VTKTOOLCURSOR_EXPORT vtkWindowLevelTool : public vtkTool
{
public:
  // Description:
  // Instantiate the object.
  static vtkWindowLevelTool *New();

  // Description:
  // Standard vtkObject methods
  vtkTypeMacro(vtkWindowLevelTool, vtkTool);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  // Description:
  // These are the methods that are called when the action takes place.
  void StartAction() override;
  void StopAction() override;
  void DoAction() override;

protected:
  vtkWindowLevelTool();
  ~vtkWindowLevelTool();

  double StartWindowLevel[2];
  vtkImageProperty *CurrentImageProperty;

  void SetCurrentImageToNthImage(int i);

private:
  vtkWindowLevelTool(const vtkWindowLevelTool&) = delete;
  void operator=(const vtkWindowLevelTool&) = delete;
};

#endif
