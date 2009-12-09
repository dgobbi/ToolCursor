/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkSurfaceCursor.h,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkSurfaceCursor - Cursor for picking and manipulating Prop3Ds
// .SECTION Description
// This class assists with picking and with providing interaction with
// objects in a 3D scene.  It allows the picking to be customized for
// different actors in the scene.

#ifndef __vtkSurfaceCursor_h
#define __vtkSurfaceCursor_h

#include "vtkObject.h"

class vtkActor;
class vtkRenderer;
class vtkMatrix4x4;
class vtkLookupTable;
class vtkDataSet;
class vtkDataSetMapper;
class vtkDataSetCollection;
class vtkPicker;
class vtkVolumePicker;
class vtkCommand;

// Cursor states, depending on what is under the cursor
#define VTK_SCURSOR_DEFAULT      0x0
#define VTK_SCURSOR_MOVEABLE     0x1
#define VTK_SCURSOR_PUSHABLE     0x2
#define VTK_SCURSOR_ROTATEABLE   0x4
#define VTK_SCURSOR_ACTOR        0x10
#define VTK_SCURSOR_VOLUME       0x20
#define VTK_SCURSOR_IMAGE_ACTOR  0x40

class VTK_EXPORT vtkSurfaceCursor : public vtkObject
{
public:
  // Description:
  // Instantiate the object.
  static vtkSurfaceCursor *New();

  // Description:
  // Standard vtkObject methods
  vtkTypeRevisionMacro(vtkSurfaceCursor,vtkObject);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Set the Renderer for this cursor to be active in.
  void SetRenderer(vtkRenderer *renderer);
  vtkRenderer *GetRenderer() { return this->Renderer; };

  // Description:
  // Set the current (x,y) display location for the cursor.
  vtkSetVector2Macro(DisplayPosition, int);
  vtkGetVector2Macro(DisplayPosition, int);

  // Description:
  // Compute the cursor position based on the display position.
  void ComputePosition();

  // Description:
  // Get the cursor position in world coordinates.
  vtkGetVector3Macro(Position, double);

  // Description:
  // Get the cursor normal in world coordinates.
  vtkGetVector3Macro(Normal, double);

  // Description:
  // Get the cursor vertical direction in world coordinates.
  // This is the ViewUp of the camera, tilted to be made perpendicular
  // to the Normal.
  vtkGetVector3Macro(Vector, double);

  // Description:
  // Get the cursor's pose matrix.  This matrix is composed from
  // the Position, Normal, and CursorVertical.
  vtkMatrix4x4 *GetMatrix() { return this->Matrix; };

  // Description:
  // Get the picker.  The picker is a vtkVolumePicker.  Every time that
  // ComputePosition is called, a pick is done and all of the picker
  // information is updated.
  vtkVolumePicker *GetPicker() { return this->Picker; };

  // Description:
  // Get the state of the cursor.  This refers to what is under the cursor,
  // and hence what kinds of actions the cursor can take.
  int GetState() { return this->State; };

  // Description:
  // Set the cursor shaped.
  void SetShape(int shape);
  int GetShape() { return this->Shape; };

  // Description:
  // Add a cursor shape.  The id for that shape will be returned.  Once
  // added, a shape cannot be removed.
  int AddShape(vtkDataSet *cursor);

  // Description:
  // Set the scale factor for the cursor, usually to make it larger.
  vtkSetMacro(Scale, double);
  vtkGetMacro(Scale, double);

  // Description:
  // Set the colors to use for the cursor.  The first color is used for
  // the part of the cursor on the "near" side, and the second color is
  // used for the part of the cursor on the "far".  Other colors are
  // specific to particular cursor shaped.  If a cursor shape uses
  // RGBA color instead of indexed color, then these colors will have
  // no effect.
  void SetColor(int i, double r, double g, double b);
  void SetColor(int i, const double rgb[3]) {
    this->SetColor(i, rgb[0], rgb[1], rgb[2]); };
  double *GetColor(int i) {
    this->GetColor(i, this->Color); return this->Color; };
  void GetColor(int i, double rgb[3]);

  // Description:
  // Set whether surface normals should always point towards the camera,
  // or whether they should point away if the cursor is on the backface
  // of the surface.  The default is "On".
  vtkSetMacro(PointNormalAtCamera, int);
  vtkBooleanMacro(PointNormalAtCamera, int);
  vtkGetMacro(PointNormalAtCamera, int);

  // Description:
  // The central event handler for the cursor.  This is primarily meant for
  // internal use.
  virtual void HandleEvent(vtkObject *object, unsigned long event, void *data);

protected:
  vtkSurfaceCursor();
  ~vtkSurfaceCursor();

  double OpacityThreshold;

  int DisplayPosition[2];
  double Position[3];
  double Normal[3];
  double Vector[3];

  double Color[3];

  int PointNormalAtCamera;
  int Shape;
  int State;
  double Scale;
  vtkMatrix4x4 *Matrix;
  vtkDataSetCollection *Shapes;
  vtkDataSetMapper *Mapper;
  vtkLookupTable *LookupTable;
  vtkActor *Actor;
  vtkVolumePicker *Picker;
  vtkRenderer *Renderer;
  vtkCommand *Command;

  virtual void ComputeState();
  virtual void MakeDefaultShapes();

  static vtkDataSet *MakePointerShape();
  static vtkDataSet *MakeCrosshairsShape();
  static vtkDataSet *MakeCrossShape(int splitCross);
  static vtkDataSet *MakeSphereShape(int splitSphere);
  static vtkDataSet *MakeConeShape(int doubleCone);

  static double ComputeScale(const double position[3], vtkRenderer *renderer);
  static void ComputeMatrix(const double position[3], const double normal[3],
                            const double vector[3], vtkMatrix4x4 *matrix); 
  static void ComputeVectorFromNormal(const double normal[3], double vector[3],
                                      vtkRenderer *renderer);

  static void UpdatePropsForPick(vtkPicker *picker, vtkRenderer *renderer);

private:
  vtkSurfaceCursor(const vtkSurfaceCursor&);  //Not implemented
  void operator=(const vtkSurfaceCursor&);  //Not implemented
};

#endif
