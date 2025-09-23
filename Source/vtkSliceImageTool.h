/*=========================================================================

  Program:   ToolCursor
  Module:    vtkSliceImageTool.h

  Copyright (c) 2010 David Gobbi
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkSliceImageTool - Move the camera focal plane in and out.
// .SECTION Description
// This class moves the focal point of the camera away from or towards the
// viewer, in order to adjust the slice plane for the images.

#ifndef vtkSliceImageTool_h
#define vtkSliceImageTool_h

#include "vtkToolCursorModule.h" // For export macro
#include "vtkImageTool.h"

class VTKTOOLCURSOR_EXPORT vtkSliceImageTool : public vtkImageTool
{
public:
  // Description:
  // Instantiate the object.
  static vtkSliceImageTool *New();

  // Description:
  // Standard vtkObject methods
  vtkTypeMacro(vtkSliceImageTool, vtkImageTool);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  // Description:
  // Turn this on to make the interaction jump to the nearest slice,
  // instead of interpolating between slices.  This option is ignored
  // if the view plane is oblique to the original image scan planes.
  // The image that is used is whatever image is on top, or if an image
  // stack is present, then the active layer of the stack is used.
  void SetJumpToNearestSlice(int val);
  void JumpToNearestSliceOff() { this->SetJumpToNearestSlice(0); }
  void JumpToNearestSliceOn() { this->SetJumpToNearestSlice(1); }
  int GetJumpToNearestSlice() { return this->JumpToNearestSlice; }

  // Description:
  // These are the methods that are called when the action takes place.
  void StartAction() override;
  void StopAction() override;
  void DoAction() override;

  // Description:
  // This is useful methods for moving back or forth by one slice.
  // The delta is the number of slices to advance by, use negative
  // values to move backwards.
  virtual void AdvanceSlice(int delta);

protected:
  vtkSliceImageTool();
  ~vtkSliceImageTool();

  int JumpToNearestSlice;
  double StartDistance;

private:
  vtkSliceImageTool(const vtkSliceImageTool&) = delete;
  void operator=(const vtkSliceImageTool&) = delete;
};

#endif
