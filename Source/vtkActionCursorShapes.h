/*=========================================================================

  Program:   ToolCursor
  Module:    vtkActionCursorShapes.h

  Copyright (c) 2010 David Gobbi
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkActionCursorShapes - Cursor shapes for interaction.
// .SECTION Description
// This class is a collection of action shapes for use with the
// vtkToolCursor.  The shapes available are "Move", "Rotate",
// "Push", and "Spin".  All cursors are constructed from the
// same 3D arrow shape and have a two-tone color scheme.

#ifndef vtkActionCursorShapes_h
#define vtkActionCursorShapes_h

#include "vtkToolCursorModule.h" // For export macro
#include "vtkCursorShapes.h"

class vtkPolyData;

class VTKTOOLCURSOR_EXPORT vtkActionCursorShapes : public vtkCursorShapes
{
public:
  // Description:
  // Instantiate the object.
  static vtkActionCursorShapes *New();

  // Description:
  // Standard vtkObject methods
  vtkTypeMacro(vtkActionCursorShapes,vtkCursorShapes);
  void PrintSelf(ostream& os, vtkIndent indent) override;

protected:
  vtkActionCursorShapes();
  ~vtkActionCursorShapes();

  void MakeShapes();
  static vtkDataSet *MakeMoveShape(vtkPolyData *arrow, int warped);
  static vtkDataSet *MakePushShape(vtkPolyData *arrow);
  static vtkDataSet *MakeSpinShape(vtkPolyData *arrow);
  static vtkDataSet *MakeZoomShape(vtkPolyData *arrow);
  static vtkPolyData *MakeArrow();
  static vtkPolyData *MakeWarpedArrow(vtkPolyData *arrow,
                                      double warpX, double warpY,
                                      double warpZ, double warpScale);

private:
  vtkActionCursorShapes(const vtkActionCursorShapes&) = delete;
  void operator=(const vtkActionCursorShapes&) = delete;
};

#endif
