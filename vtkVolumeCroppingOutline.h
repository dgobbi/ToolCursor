/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkVolumeCroppingOutline.h,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkVolumeCroppingOutline - wireframe outline for volume cropping
// .SECTION Description
// vtkVolumeCroppingOutline generates a wireframe outline that corresponds
// to the cropping region of a vtkVolumeMapper.  

#ifndef __vtkVolumeCroppingOutline_h
#define __vtkVolumeCroppingOutline_h

#include "vtkPolyDataAlgorithm.h"

class vtkVolumeMapper;

class VTK_EXPORT vtkVolumeCroppingOutline : public vtkPolyDataAlgorithm
{
public:
  static vtkVolumeCroppingOutline *New();
  vtkTypeRevisionMacro(vtkVolumeCroppingOutline,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Set the mapper that has the cropping region that the outline will
  // be generated for.  The mapper must have an input, because the
  // bounds of the data must be computed in order to generate the
  // outline.
  virtual void SetVolumeMapper(vtkVolumeMapper *mapper);
  vtkVolumeMapper *GetVolumeMapper() { return this->VolumeMapper; };

  // Description:
  // Set whether to add color scalars to the lines.  By default,
  // the output has no scalars and the color must be set in the
  // property of the actor.
  vtkSetMacro(UseColorScalars, int);
  vtkBooleanMacro(UseColorScalars, int);
  vtkGetMacro(UseColorScalars, int);

  // Description:
  // Set the color of the outline.  This has no effect unless ColorOutline
  // is On.  The default color is red.
  vtkSetVector3Macro(Color, double);
  vtkGetVector3Macro(Color, double);

  // Description:
  // Set the active plane, e.g. to display which plane is currently being
  // modified by an interaction.  Set this to -1 if there is no active plane.
  // The default value is -1.
  vtkSetMacro(ActivePlaneId, int);
  vtkGetMacro(ActivePlaneId, int);

  // Description:
  // Set the color of the active cropping plane.  This has no effect unless
  // ColorOutline is On and ActivePlaneId is non-negative.  Default is yellow. 
  vtkSetVector3Macro(ActivePlaneColor, double);
  vtkGetVector3Macro(ActivePlaneColor, double);

protected:
  vtkVolumeCroppingOutline();
  ~vtkVolumeCroppingOutline();

  vtkVolumeMapper *VolumeMapper;
  int UseColorScalars;
  int ActivePlaneId;
  double Color[3];
  double ActivePlaneColor[3];

  int Cropping;
  int CroppingRegionFlags;
  double Bounds[6];
  double CroppingRegionPlanes[6];

  int ComputeCubePlanes(double planes[3][4]);

  virtual int ComputePipelineMTime(vtkInformation* request,
                                   vtkInformationVector** inputVector,
                                   vtkInformationVector* outputVector,
                                   int requestFromOutputPort,
                                   unsigned long* mtime);

  virtual int RequestInformation(vtkInformation* request,
                                 vtkInformationVector** inputVector,
                                 vtkInformationVector* outputVector);

  virtual int RequestData(vtkInformation* request,
                          vtkInformationVector** inputVector,
                          vtkInformationVector* outputVector);

private:
  vtkVolumeCroppingOutline(const vtkVolumeCroppingOutline&);  // Not implemented.
  void operator=(const vtkVolumeCroppingOutline&);  // Not implemented.
};

#endif
