/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkSurfaceCursorAction.h,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkSurfaceCursorAction - A base class for an interaction state.
// .SECTION Description
// This is the base class for interactions, i.e. for what happens when
// you drag the mouse in a renderer.  It contains all the necessary
// state information, data, and code for the interaction.

#ifndef __vtkSurfaceCursorAction_h
#define __vtkSurfaceCursorAction_h

#include "vtkObject.h"

class vtkSurfaceCursor;

class VTK_EXPORT vtkSurfaceCursorAction : public vtkObject
{
public:
  // Description:
  // Instantiate the object.
  static vtkSurfaceCursorAction *New();

  // Description:
  // Standard vtkObject methods
  vtkTypeRevisionMacro(vtkSurfaceCursorAction,vtkObject);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Set the cursor that this action is bound to.
  void SetSurfaceCursor(vtkSurfaceCursor *cursor);
  vtkSurfaceCursor *GetSurfaceCursor() { return this->SurfaceCursor; };

  // Description:
  // These are the methods that are called when the action takes place.
  virtual void StartAction();
  virtual void StopAction();
  virtual void DoAction();

  // Description:
  // This method allow the action to constrain the cursor's position.
  // It is called immediately after the cursor positioning pick occurs.
  // The default action is to do nothing.
  virtual void ConstrainCursor(double position[3], double normal[3]);

  // Description:
  // Get the current display position.
  void GetDisplayPosition(double &x, double &y) {
    x = this->DisplayPosition[0]; y = this->DisplayPosition[1]; };
  double *GetDisplayPosition() { return this->DisplayPosition; };

  // Description:
  // Get the previous display position.
  void GetLastDisplayPosition(double &x, double &y) {
    x = this->LastDisplayPosition[0]; y = this->LastDisplayPosition[1]; };
  double *GetLastDisplayPosition() { return this->LastDisplayPosition; };

  // Description:
  // Get the display position from the start of the action.
  void GetStartDisplayPosition(double &x, double &y) {
    x = this->StartDisplayPosition[0]; y = this->StartDisplayPosition[1]; };
  double *GetStartDisplayPosition() { return this->StartDisplayPosition; };

  // Description:
  // Get the world coordinates of the cursor at the start of the action.
  void GetStartPosition(double position[3]) {
    position[0] = this->StartPosition[0];
    position[1] = this->StartPosition[1];
    position[2] = this->StartPosition[2]; };
  double *GetStartPosition() { return this->StartPosition; };

protected:
  vtkSurfaceCursorAction();
  ~vtkSurfaceCursorAction();

  vtkSurfaceCursor *SurfaceCursor;

  double StartPosition[3];
  double StartDisplayPosition[2];
  double LastDisplayPosition[2];
  double DisplayPosition[2];

  // Description:
  // Convert from world coords to display coords.
  void WorldToDisplay(const double world[3], double &x, double &y, double &z);

  // Description:
  // Convert from display coords to world coords.
  void DisplayToWorld(double x, double y, double z, double world[3]);

  // Description:
  // Given an (x,y,z) display coord, provide the corresponding world-space
  // point and the vector along the view ray for that point.
  void GetViewRay(double x, double y, double z, double p[3], double v[3]);

private:
  vtkSurfaceCursorAction(const vtkSurfaceCursorAction&);  //Not implemented
  void operator=(const vtkSurfaceCursorAction&);  //Not implemented
};

#endif
