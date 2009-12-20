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

protected:
  vtkSurfaceCursorAction();
  ~vtkSurfaceCursorAction();

  vtkSurfaceCursor *SurfaceCursor;

  double StartPosition[3];
  double StartDisplayPosition[2];
  double LastDisplayPosition[2];
  double DisplayPosition[2];

private:
  vtkSurfaceCursorAction(const vtkSurfaceCursorAction&);  //Not implemented
  void operator=(const vtkSurfaceCursorAction&);  //Not implemented
};

#endif
