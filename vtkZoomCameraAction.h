/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkZoomCameraAction.h,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkZoomCameraAction - Zoom the view by moving the camera in and out.
// .SECTION Description
// This class is used to zoom the view.  The default way of doing the zoom
// is by dollying the camera towards or away from the focal point in order
// to achieve the zoom.

#ifndef __vtkZoomCameraAction_h
#define __vtkZoomCameraAction_h

#include "vtkSurfaceCursorAction.h"

class vtkTransform;

class VTK_EXPORT vtkZoomCameraAction : public vtkSurfaceCursorAction
{
public:
  // Description:
  // Instantiate the object.
  static vtkZoomCameraAction *New();

  // Description:
  // Standard vtkObject methods
  vtkTypeRevisionMacro(vtkZoomCameraAction,vtkSurfaceCursorAction);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Set whether to zoom by dollying the camera or by decreasing the
  // field of view.  This option is ignored if the camera is using
  // a parallel projection.  The default is to zoom by dollying.
  vtkSetMacro(ZoomByDolly, int);
  vtkBooleanMacro(ZoomByDolly, int);
  int GetZoomByDolly(int) { return this->ZoomByDolly; };

  // Description:
  // These are the methods that are called when the action takes place.
  virtual void StartAction();
  virtual void StopAction();
  virtual void DoAction();

  // Description:
  // Constrain the position or orientation of the cursor.
  virtual void ConstrainCursor(double position[3], double normal[3]);

protected:
  vtkZoomCameraAction();
  ~vtkZoomCameraAction();

  int ZoomByDolly;

  double StartCameraPosition[3];
  double StartClippingRange[2];
  double StartParallelScale;
  double StartViewAngle;

  double ZoomFactor;

  vtkTransform *Transform;

private:
  vtkZoomCameraAction(const vtkZoomCameraAction&);  //Not implemented
  void operator=(const vtkZoomCameraAction&);  //Not implemented
};

#endif
