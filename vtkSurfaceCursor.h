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
class vtkRenderWindowInteractor;
class vtkMatrix4x4;
class vtkLookupTable;
class vtkPolyData;
class vtkDataSet;
class vtkDataSetMapper;
class vtkDataSetCollection;
class vtkCollection;
class vtkPicker;
class vtkVolumePicker;
class vtkCommand;

class vtkSurfaceCursorAction;

// Cursor actions.
#define VTK_SCURSOR_PUSH          1
#define VTK_SCURSOR_ROTATE        2
#define VTK_SCURSOR_SPIN          3

// Event modifiers, which usually control the mode.
#define VTK_SCURSOR_SHIFT        0x01
#define VTK_SCURSOR_CAPS         0x02
#define VTK_SCURSOR_CONTROL      0x04
#define VTK_SCURSOR_META         0x08
#define VTK_SCURSOR_ALT          0x16
#define VTK_SCURSOR_B1          0x100
#define VTK_SCURSOR_B2          0x200
#define VTK_SCURSOR_B3          0x400

// Pick flags, these describe what is under the cursor.
#define VTK_SCURSOR_PROP3D       0x0F00
#define VTK_SCURSOR_ACTOR        0x0100
#define VTK_SCURSOR_VOLUME       0x0200
#define VTK_SCURSOR_IMAGE_ACTOR  0x0400
#define VTK_SCURSOR_CLIP_PLANE   0x1000
#define VTK_SCURSOR_CROP_PLANE   0x2000

// Cursor shapes. Copies of basic system cursors
#define VTK_SCURSOR_POINTER      0
#define VTK_SCURSOR_CROSSHAIRS   1
// Cursor shapes. Geometrical shapes.
#define VTK_SCURSOR_CROSS        2
#define VTK_SCURSOR_CROSS_SPLIT  3
#define VTK_SCURSOR_CONE         4
#define VTK_SCURSOR_CONE_DUAL    5
#define VTK_SCURSOR_SPHERE       6
#define VTK_SCURSOR_SPHERE_SPLIT 7
// Cursor shapes. Action cursors.
#define VTK_SCURSOR_MOVER        8
#define VTK_SCURSOR_ROCKER       9
#define VTK_SCURSOR_PUSHER       10
#define VTK_SCURSOR_SPINNER      11

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
  vtkSetVector2Macro(DisplayPosition, double);
  vtkGetVector2Macro(DisplayPosition, double);

  // Description:
  // Compute the cursor position based on the display position.
  virtual void ComputePosition();

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
  // Get the pick flags.  The flags provide information about what was
  // under the cursor the last time that a pick was done.  The flags are
  // used to determine what kinds of actions the cursor can take.
  int GetPickFlags() { return this->PickFlags; };

  // Description:
  // Set the current mode for the cursor.  The mode is usually computed
  // from the Modifier ivar by the ComputeMode() method, i.e. the mode will
  // usually change in response to modifier keys like "Shift" and "Control".
  // If you don't use SetModifier() or BindInteractor(), then you can set
  // the mode manually with this method.
  void SetMode(int mode);
  int GetMode() { return this->Mode; };

  // Description:
  // Set the cursor shaped.  If you use BindInteractor() or SetModifier()
  // then the shape will be computed from the Mode and the PickFlags.  If
  // you don't use these methods, then you can set the shape manually.
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
  // Bind an interactor.  This will add observers for all mouse and keyboard
  // events that the cursor needs.  All bound events are observed via the
  // HandleEvent() method.  An alternative to binding an Interactor is to
  // call the MoveToDisplayPosition(), SetModifier(), and SetMouseInRenderer()
  // methods directly.
  virtual void BindInteractor(vtkRenderWindowInteractor *iren);

  // Description:
  // The central event handler for the cursor.  All interactor events go
  // through this method. 
  virtual void HandleEvent(vtkObject *object, unsigned long event, void *data);

  // Description:
  // Move the cursor to a specific position, usually in response to the
  // mouse motion.  This is not a passive method like SetDisplayPosition().
  // Depending on the Action ivar, the motion will go to the appropriate
  // interaction method. 
  virtual void MoveToDisplayPosition(double x, double y);

  // Description:
  // Set or get the modifier bitfield.  This will set the Mode and the
  // cursor Shape, and if the modifier bits indicate that a mouse button
  // has been pressed, it will also set the Action.  If you want to set
  // the Action manually (i.e. if you don't like the default bindings)
  // then you should not call this method.
  virtual void SetModifier(int modifier);  
  int GetModifier() { return this->Modifier; };

  // Description:
  // Set the current action.  If you use BindInteractor() or SetModifier(),
  // then this will be set automatically in response to button events.  If
  // you don't use either of these methods, then you can set Action
  // manually.  The value of Action controls what happens when
  // MoveToDisplayPosition() is called.
  virtual void SetAction(int action);
  int GetAction() { return this->Action; };

  // Description:
  // Add an action.  The id for the action will be returned.  Once an
  // action is added, it cannot be removed.
  int AddAction(vtkSurfaceCursorAction *action);

  // Description:
  // Set whether the mouse is in the renderer.  This controls cursor
  // visibility.  If you use BindToInteractor(), this method is called
  // automatically when required.
  virtual void SetMouseInRenderer(int inside);
  int GetMouseInRenderer() { return this->MouseInRenderer; };

  // Description:
  // We override this method to modify the actor, otherwise the
  // RenderWindow won't know that it needs to render.
  virtual void Modified();

protected:
  vtkSurfaceCursor();
  ~vtkSurfaceCursor();

  double OpacityThreshold;

  double DisplayPosition[2];
  double Position[3];
  double Normal[3];
  double Vector[3];

  double Color[3];

  int PointNormalAtCamera;
  int Shape;
  int Mode;
  int PickFlags;
  int Action;
  int ActionButton;
  int Modifier;
  int MouseInRenderer;
  double Scale;
  vtkMatrix4x4 *Matrix;
  vtkDataSetCollection *Shapes;
  vtkCollection *Actions;
  vtkDataSetMapper *Mapper;
  vtkLookupTable *LookupTable;
  vtkActor *Actor;
  vtkVolumePicker *Picker;
  vtkRenderer *Renderer;
  vtkCommand *Command;

  virtual int ComputeMode(int modifier);
  virtual int ComputeShape(int mode, int pickFlags);
  virtual int ComputeAction(int mode, int pickFlags, int button);

  virtual void MakeDefaultShapes();
  static vtkDataSet *MakePointerShape();
  static vtkDataSet *MakeCrosshairsShape();
  static vtkDataSet *MakeCrossShape(int splitCross);
  static vtkDataSet *MakeSphereShape(int splitSphere);
  static vtkDataSet *MakeConeShape(int doubleCone);
  static vtkDataSet *MakeMoverShape(int warped);
  static vtkDataSet *MakePusherShape();
  static vtkDataSet *MakeSpinnerShape();
  static vtkPolyData *MakeWarpedArrow(double warpX, double warpY,
                                      double warpZ, double warpScale);

  static int ComputePickFlags(vtkVolumePicker *picker);
  static double ComputeScale(const double position[3], vtkRenderer *renderer);
  static void ComputeMatrix(const double position[3], const double normal[3],
                            const double vector[3], vtkMatrix4x4 *matrix); 
  static void ComputeVectorFromNormal(const double normal[3], double vector[3],
                                      vtkDataSetMapper *cursorMapper,
                                      vtkRenderer *renderer);

  static void UpdatePropsForPick(vtkPicker *picker, vtkRenderer *renderer);

  static int ModifierFromKeySym(const char *keysym);

private:
  vtkSurfaceCursor(const vtkSurfaceCursor&);  //Not implemented
  void operator=(const vtkSurfaceCursor&);  //Not implemented
};

#endif
