/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkSurfaceCursorInteractorObserver.h,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkSurfaceCursorInteractorObserver - Observer for vtkSurfaceCursor
// .SECTION Description
// This is a helper class for vtkSurfaceCursor.  Specifically, this class
// receives events from a vtkRenderWindowInteractor and then calls the
// appropriate methods of vtkSurfaceCursor.

#ifndef __vtkSurfaceCursorInteractorObserver_h
#define __vtkSurfaceCursorInteractorObserver_h

#include "vtkInteractorObserver.h"

class vtkSurfaceCursor;

class VTK_EXPORT vtkSurfaceCursorInteractorObserver :
  public vtkInteractorObserver
{
public:
  // Description:
  // Instantiate the object.
  static vtkSurfaceCursorInteractorObserver *New();

  // Description:
  // Standard vtkObject methods
  vtkTypeRevisionMacro(vtkSurfaceCursorInteractorObserver,
                       vtkInteractorObserver);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Enable the cursor.
  virtual void SetEnabled(int enable);

  // Description:
  // Set the cursor that this object observes events for.
  void SetSurfaceCursor(vtkSurfaceCursor *cursor) {
     this->SurfaceCursor = cursor; };
  vtkSurfaceCursor *GetSurfaceCursor() {
     return this->SurfaceCursor; };

  // Description:
  // Get vtkSurfaceCursor "modifier" bits from a VTK keysym.
  static int ModifierFromKeySym(const char *keysym);

protected:
  vtkSurfaceCursorInteractorObserver();
  ~vtkSurfaceCursorInteractorObserver();

  static void ProcessEvents(vtkObject* object,
                            unsigned long event,
                            void* clientdata,
                            void* calldata);

  vtkSurfaceCursor *SurfaceCursor;

private:
  vtkSurfaceCursorInteractorObserver(const vtkSurfaceCursorInteractorObserver&);  //Not implemented
  void operator=(const vtkSurfaceCursorInteractorObserver&);  //Not implemented
};

#endif
