/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkPanCameraAction.h,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkPanCameraAction - Pan the camera in order to move the scene.
// .SECTION Description
// This class is used to pan the camera.  The pan is achieved by moving
// the camera position and focal point together, not by slowly rotating
// the camera around its position.

#ifndef __vtkPanCameraAction_h
#define __vtkPanCameraAction_h

#include "vtkSurfaceCursorAction.h"

class vtkTransform;

class VTK_EXPORT vtkPanCameraAction : public vtkSurfaceCursorAction
{
public:
  // Description:
  // Instantiate the object.
  static vtkPanCameraAction *New();

  // Description:
  // Standard vtkObject methods
  vtkTypeRevisionMacro(vtkPanCameraAction,vtkSurfaceCursorAction);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // These are the methods that are called when the action takes place.
  virtual void StartAction();
  virtual void StopAction();
  virtual void DoAction();

protected:
  vtkPanCameraAction();
  ~vtkPanCameraAction();

  double StartCameraFocalPoint[3];
  double StartCameraPosition[3];

  vtkTransform *Transform;

private:
  vtkPanCameraAction(const vtkPanCameraAction&);  //Not implemented
  void operator=(const vtkPanCameraAction&);  //Not implemented
};

#endif
