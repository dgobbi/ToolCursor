/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkRotateCameraAction.h,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkRotateCameraAction - Controls camera rotation.
// .SECTION Description
// This class controls camera rotation interaction.

#ifndef __vtkRotateCameraAction_h
#define __vtkRotateCameraAction_h

#include "vtkSurfaceCursorAction.h"

class vtkTransform;

class VTK_EXPORT vtkRotateCameraAction : public vtkSurfaceCursorAction
{
public:
  // Description:
  // Instantiate the object.
  static vtkRotateCameraAction *New();

  // Description:
  // Standard vtkObject methods
  vtkTypeRevisionMacro(vtkRotateCameraAction,vtkSurfaceCursorAction);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // These are the methods that are called when the action takes place.
  virtual void StartAction();
  virtual void StopAction();
  virtual void DoAction();

  // Description:
  // This method allows the action to constrain the cursor position.
  virtual void ConstrainCursor(double position[3], double normal[3]);

protected:
  vtkRotateCameraAction();
  ~vtkRotateCameraAction();

  void SetNormal(double normal[3]) {
    this->Normal[0] = normal[0];
    this->Normal[1] = normal[1];
    this->Normal[2] = normal[2]; };

  void GetNormal(double normal[3]) {
    normal[0] = this->Normal[0];
    normal[1] = this->Normal[1];
    normal[2] = this->Normal[2]; };

  void GetCenterOfRotation(double center[3]) {
    center[0] = this->CenterOfRotation[0];
    center[1] = this->CenterOfRotation[1];
    center[2] = this->CenterOfRotation[2]; };

  double GetRadius() {
    return this->Radius; };

  double Normal[3];
  double CenterOfRotation[3];
  double Radius;

  double StartCameraPosition[3];
  double StartCameraViewUp[3];

  vtkTransform *Transform;

private:
  vtkRotateCameraAction(const vtkRotateCameraAction&);  //Not implemented
  void operator=(const vtkRotateCameraAction&);  //Not implemented
};

#endif
