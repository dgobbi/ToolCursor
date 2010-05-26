/*=========================================================================

  Program:   ToolCursor
  Module:    vtkVolumeCroppingOutline.h

  Copyright (c) 2010 David Gobbi
  All rights reserved.

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
  vtkTypeMacro(vtkVolumeCroppingOutline,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Set the mapper that has the cropping region that the outline will
  // be generated for.  The mapper must have an input, because the
  // bounds of the data must be computed in order to generate the
  // outline.
  virtual void SetVolumeMapper(vtkVolumeMapper *mapper);
  vtkVolumeMapper *GetVolumeMapper() { return this->VolumeMapper; };

  // Description:
  // Set whether to add color scalars to the output.  By default,
  // the output has no scalars and the color must be set in the
  // property of the actor.
  vtkSetMacro(GenerateScalars, int);
  vtkBooleanMacro(GenerateScalars, int);
  vtkGetMacro(GenerateScalars, int);

  // Description:
  // Set whether to generate an outline wherever an input face was
  // cut by a plane.  This is on by default.
  vtkSetMacro(GenerateOutline, int);
  vtkBooleanMacro(GenerateOutline, int);
  vtkGetMacro(GenerateOutline, int);

  // Description:
  // Set whether to generate polygonal faces for the output.  By default,
  // only lines are generated.
  vtkSetMacro(GenerateFaces, int);
  vtkBooleanMacro(GenerateFaces, int);
  vtkGetMacro(GenerateFaces, int);

  // Description:
  // Set the color of the outline.  This has no effect unless GenerateScalars
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
  int GenerateScalars;
  int GenerateOutline;
  int GenerateFaces;
  int ActivePlaneId;
  double Color[3];
  double ActivePlaneColor[3];

  int Cropping;
  int CroppingRegionFlags;
  double Bounds[6];
  double CroppingRegionPlanes[6];

  static int ComputeCubePlanes(double planes[3][4],
                               double croppingPlanes[6],
                               double bounds[6]);

  static void GeneratePolys(vtkCellArray *polys,
                            vtkUnsignedCharArray *scalars,
                            unsigned char colors[2][3],
                            int activePlane,
                            int flags,
                            int tolPtId[3][4]);

  static void GenerateLines(vtkCellArray *lines,
                            vtkUnsignedCharArray *scalars,
                            unsigned char colors[2][3],
                            int activePlane,
                            int flags,
                            int tolPtId[3][4]);

  static void GeneratePoints(vtkPoints *points,
                             vtkCellArray *lines,
                             vtkCellArray *polys,
                             double planes[3][4],
                             double tol);

  static void NudgeCropPlanesToBounds(int tolPtId[3][4],
                                      double planes[3][4],
                                      double tol);

  static void CreateColorValues(unsigned char colors[2][3],
                                double color1[3], double color2[3]);

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
