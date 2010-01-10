/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkSurfaceCursorShapes.h,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkSurfaceCursorShapes - A base class for cursor shape sets.
// .SECTION Description
// This is the base class for vtkSurfaceCursor cursor shapes.  Each
// subclass will define a set of 3D cursor shapes that can be used.

#ifndef __vtkSurfaceCursorShapes_h
#define __vtkSurfaceCursorShapes_h

#include "vtkObject.h"

class vtkSurfaceCursorShapeArray;
class vtkDataSet;

// Flags for shapes.
#define VTK_SCURSOR_FLATX    0x01  // cursor is mainly in YZ plane
#define VTK_SCURSOR_FLATY    0x02  // cursor is mainly in XZ plane
#define VTK_SCURSOR_RGB      0x10  // cursor uses RGB scalars
#define VTK_SCURSOR_NORMLOCK 0x40  // lock orientation during interaction

class VTK_EXPORT vtkSurfaceCursorShapes : public vtkObject
{
public:
  // Description:
  // Instantiate the object.
  static vtkSurfaceCursorShapes *New();

  // Description:
  // Standard vtkObject methods
  vtkTypeRevisionMacro(vtkSurfaceCursorShapes,vtkObject);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Get the number of shapes.
  int GetNumberOfShapes() { return this->NumberOfShapes; };

  // Description:
  // Add a cursor shape.  The id for that shape will be returned.  Once
  // added, a shape cannot be changed or removed.
  int AddShape(const char *name, vtkDataSet *shape, int flags);

  // Description:
  // Get the name of the shape at the specified index.  Returns null if
  // the index is out of range.
  const char *GetShapeName(int i);

  // Description:
  // Get the index of the shape with the specified name.  Returns -1 if the
  // name is not recognized.
  int GetShapeIndex(const char *name);

  // Description:
  // Get the specified shape.  Returns null if the index is out of range or
  // or if the name is not recognized.
  vtkDataSet *GetShapeData(int i);
  vtkDataSet *GetShapeData(const char *name) {
    return this->GetShapeData(this->GetShapeIndex(name)); };

  // Description:
  // Get the flags for this shape.  Returns zero if the index is out of range
  // or if the name is not recognized.
  int GetShapeFlags(int i);
  int GetShapeFlags(const char *name) {
    return this->GetShapeFlags(this->GetShapeIndex(name)); };

protected:
  vtkSurfaceCursorShapes();
  ~vtkSurfaceCursorShapes();

  int NumberOfShapes;
  vtkSurfaceCursorShapeArray *Shapes;

private:
  vtkSurfaceCursorShapes(const vtkSurfaceCursorShapes&);  //Not implemented
  void operator=(const vtkSurfaceCursorShapes&);  //Not implemented
};

#endif
