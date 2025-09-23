/*=========================================================================

  Program:   ToolCursor
  Module:    vtkSpinCameraTool.h

  Copyright (c) 2010 David Gobbi
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkSpinCameraTool - Spin the camera with the cursor.
// .SECTION Description
// This is an interaction class for spinning the cursor around the
// camera's view direction.

#ifndef vtkSpinCameraTool_h
#define vtkSpinCameraTool_h

#include "vtkToolCursorModule.h" // For export macro
#include "vtkTool.h"

class vtkTransform;

class VTKTOOLCURSOR_EXPORT vtkSpinCameraTool : public vtkTool
{
public:
  // Description:
  // Instantiate the object.
  static vtkSpinCameraTool *New();

  // Description:
  // Standard vtkObject methods
  vtkTypeMacro(vtkSpinCameraTool,vtkTool);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  // Description:
  // These are the methods that are called when the action takes place.
  void StartAction() override;
  void StopAction() override;
  void DoAction() override;

protected:
  vtkSpinCameraTool();
  ~vtkSpinCameraTool();

  double StartCameraViewUp[3];

  vtkTransform *Transform;

private:
  vtkSpinCameraTool(const vtkSpinCameraTool&) = delete;
  void operator=(const vtkSpinCameraTool&) = delete;
};

#endif
