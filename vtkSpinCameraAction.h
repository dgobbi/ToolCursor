/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkSpinCameraAction.h,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkSpinCameraAction - Spin the camera with the cursor.
// .SECTION Description
// This is an interaction class for spinning the cursor around the
// camera's view direction.

#ifndef __vtkSpinCameraAction_h
#define __vtkSpinCameraAction_h

#include "vtkSurfaceCursorAction.h"

class vtkTransform;

class VTK_EXPORT vtkSpinCameraAction : public vtkSurfaceCursorAction
{
public:
  // Description:
  // Instantiate the object.
  static vtkSpinCameraAction *New();

  // Description:
  // Standard vtkObject methods
  vtkTypeRevisionMacro(vtkSpinCameraAction,vtkSurfaceCursorAction);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // These are the methods that are called when the action takes place.
  virtual void StartAction();
  virtual void StopAction();
  virtual void DoAction();

protected:
  vtkSpinCameraAction();
  ~vtkSpinCameraAction();

  double StartCameraViewUp[3];

  vtkTransform *Transform;

private:
  vtkSpinCameraAction(const vtkSpinCameraAction&);  //Not implemented
  void operator=(const vtkSpinCameraAction&);  //Not implemented
};

#endif
