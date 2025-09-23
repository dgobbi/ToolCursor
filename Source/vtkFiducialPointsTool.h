/*=========================================================================

  Program:   ToolCursor
  Module:    vtkFiducialPointsTool.h

  Copyright (c) 2010 David Gobbi
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkFiducialPointsTool - Display and control fiducial markers
// .SECTION Description
// This class allows the placement of fiducial markers.

#ifndef vtkFiducialPointsTool_h
#define vtkFiducialPointsTool_h

#include "vtkTool.h"

class vtkTransform;
class vtkGlyph3D;
class vtkPoints;
class vtkPolyData;
class vtkPointSet;
class vtkDataSetMapper;
class vtkActor;

class VTK_EXPORT vtkFiducialPointsTool : public vtkTool
{
public:
  // Description:
  // Instantiate the object.
  static vtkFiducialPointsTool *New();

  // Description:
  // Standard vtkObject methods
  vtkTypeMacro(vtkFiducialPointsTool,vtkTool);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  // Description:
  // Set the points that will used for the fiducial positions.
  virtual void SetPoints(vtkPoints *points);
  virtual vtkPoints *GetPoints();

  // Description:
  // Set the marker to use at each point.
  virtual void SetMarker(vtkPolyData *data);
  virtual vtkPolyData *GetMarker();

  // Description:
  // Get the actor for the fiducials.
  vtkActor *GetActor() { return this->Actor; };

  // Description:
  // These are the methods that are called when the action takes place.
  void StartAction() override;
  void StopAction() override;
  void DoAction() override;

protected:
  vtkFiducialPointsTool();
  ~vtkFiducialPointsTool();

  vtkPolyData *PointSet;
  vtkGlyph3D *Glyph3D;
  vtkDataSetMapper *Mapper;
  vtkActor *Actor;
  vtkTransform *Transform;

private:
  vtkFiducialPointsTool(const vtkFiducialPointsTool&) = delete;
  void operator=(const vtkFiducialPointsTool&) = delete;
};

#endif
