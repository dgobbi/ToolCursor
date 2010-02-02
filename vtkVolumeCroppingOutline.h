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

protected:
  vtkVolumeCroppingOutline();
  ~vtkVolumeCroppingOutline();

  vtkVolumeMapper *VolumeMapper;
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
