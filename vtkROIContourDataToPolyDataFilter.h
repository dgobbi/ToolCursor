/*=========================================================================

  Program:   ToolCursor
  Module:    vtkROIContourDataToPolyDataFilter.h

  Copyright (c) 2010 David Gobbi
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkROIContourDataToPolyDataFilter - Convert ROI contours to polydata
// .SECTION Description
// This filter will convert a set of contours that define an ROI into a
// vtkPolyData consisting of verts and lines.

#ifndef __vtkROIContourDataToPolyDataFilter_h
#define __vtkROIContourDataToPolyDataFilter_h

#include "vtkPolyDataAlgorithm.h"

class vtkROIContourData;
class vtkPlane;

class VTK_EXPORT vtkROIContourDataToPolyDataFilter : public vtkPolyDataAlgorithm
{
public:
  static vtkROIContourDataToPolyDataFilter *New();
  vtkTypeMacro(vtkROIContourDataToPolyDataFilter,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Set the plane to use for selecting the contours.  Only contours within
  // a certain distance of this plane will be included in the output.  If the
  // plane is set to NULL, then the output will contain all the contours.
  void SetSelectionPlane(vtkPlane *plane);
  vtkPlane *GetSelectionPlane() { return this->SelectionPlane; }

  // Description:
  // Set the tolerance for how far the contours can be from the selection plane
  // and still be included in the output.
  vtkSetMacro(SelectionPlaneTolerance, double);
  vtkGetMacro(SelectionPlaneTolerance, double);

protected:
  vtkROIContourDataToPolyDataFilter();
  ~vtkROIContourDataToPolyDataFilter();

  virtual int ComputePipelineMTime(
    vtkInformation* request, vtkInformationVector** inputVector,
    vtkInformationVector* outputVector, int requestFromOutputPort,
    unsigned long* mtime);

  virtual int RequestData(
    vtkInformation* request, vtkInformationVector** inputVector,
    vtkInformationVector* outputVector);

  virtual int FillInputPortInformation(int port, vtkInformation *info);

  vtkPlane *SelectionPlane;
  double SelectionPlaneTolerance;

private:
  vtkROIContourDataToPolyDataFilter(const vtkROIContourDataToPolyDataFilter&);  // Not implemented.
  void operator=(const vtkROIContourDataToPolyDataFilter&);  // Not implemented.
};

#endif
