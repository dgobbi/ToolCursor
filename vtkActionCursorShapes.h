/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkActionCursorShapes.h,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkActionCursorShapes - Cursor shapes for interaction.
// .SECTION Description
// This class is a collection of action shapes for use with the
// vtkSurfaceCursor.  The shapes available are "Mover", "Rotator",
// "Pusher", and "Spinner".  All cursors are constructed from the
// same 3D arrow shape and have a two-tone color scheme.

#ifndef __vtkActionCursorShapes_h
#define __vtkActionCursorShapes_h

#include "vtkSurfaceCursorShapes.h"

class vtkPolyData;

class VTK_EXPORT vtkActionCursorShapes : public vtkSurfaceCursorShapes
{
public:
  // Description:
  // Instantiate the object.
  static vtkActionCursorShapes *New();

  // Description:
  // Standard vtkObject methods
  vtkTypeRevisionMacro(vtkActionCursorShapes,vtkSurfaceCursorShapes);
  void PrintSelf(ostream& os, vtkIndent indent);

protected:
  vtkActionCursorShapes();
  ~vtkActionCursorShapes();

  void MakeShapes();
  static vtkDataSet *MakeMoverShape(vtkPolyData *arrow, int warped);
  static vtkDataSet *MakePusherShape(vtkPolyData *arrow);
  static vtkDataSet *MakeSpinnerShape(vtkPolyData *arrow);
  static vtkPolyData *MakeArrow();
  static vtkPolyData *MakeWarpedArrow(vtkPolyData *arrow,
                                      double warpX, double warpY,
                                      double warpZ, double warpScale);

private:
  vtkActionCursorShapes(const vtkActionCursorShapes&);  //Not implemented
  void operator=(const vtkActionCursorShapes&);  //Not implemented
};

#endif
