/*=========================================================================

  Program:   ToolCursor
  Module:    vtkLassoImageTool.h

  Copyright (c) 2010 David Gobbi
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkLassoImageTool - Draw contours on an image.
// .SECTION Description
// This class allows contours to be drawn on an image for manual segmentation
// and region-of-interest definition.

#ifndef __vtkLassoImageTool_h
#define __vtkLassoImageTool_h

#include "vtkImageTool.h"

class vtkMatrix4x4;
class vtkROIContourData;
class vtkROIContourDataToPolyData;
class vtkGlyph3D;
class vtkPoints;
class vtkPolyData;
class vtkPointSet;
class vtkDataSetMapper;
class vtkCellLocator;
class vtkActor;
class vtkRenderer;

class VTK_EXPORT vtkLassoImageTool : public vtkImageTool
{
public:
  // Description:
  // Instantiate the object.
  static vtkLassoImageTool *New();

  // Description:
  // Standard vtkObject methods
  vtkTypeMacro(vtkLassoImageTool,vtkImageTool);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // The ROI data to interact with.
  void SetROIContourData(vtkROIContourData *data);
  vtkROIContourData *GetROIContourData() {
    return this->ROIData; }

  // Description:
  // Set the marker to use at each point.
  virtual void SetMarker(vtkPolyData *data);
  virtual vtkPolyData *GetMarker();

  // Description:
  // These are the methods that are called when the action takes place.
  virtual void StartAction();
  virtual void StopAction();
  virtual void DoAction();

  // Description:
  // Add any associated items to the renderer.
  virtual void AddViewPropsToRenderer(vtkRenderer *renderer);
  virtual void RemoveViewPropsFromRenderer(vtkRenderer *renderer);

protected:
  vtkLassoImageTool();
  ~vtkLassoImageTool();

  vtkCellLocator *CellLocator;
  vtkROIContourData *ROIData;
  vtkROIContourDataToPolyData *ROIDataToPointSet;
  vtkROIContourDataToPolyData *ROIDataToPolyData;
  vtkGlyph3D *Glyph3D;
  vtkDataSetMapper *GlyphMapper;
  vtkActor *GlyphActor;
  vtkDataSetMapper *ContourMapper;
  vtkActor *ContourActor;
  vtkMatrix4x4 *Matrix;

  vtkIdType CurrentPointId;
  int CurrentContourId;
  double InitialPointPosition[3];

private:
  vtkLassoImageTool(const vtkLassoImageTool&);  //Not implemented
  void operator=(const vtkLassoImageTool&);  //Not implemented
};

#endif
