/*=========================================================================

  Program:   ToolCursor
  Module:    vtkClipOutlineWithPlanes.h

  Copyright (c) 2010 David Gobbi
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkClipOutlineWithPlanes - Clip an outline with a plane collection
// .SECTION Description
// vtkClipOutlineWithPlanes will clip a closed polydata surface with a
// collection of clipping planes.  Its default action is to generate
// an outline wherever the clipping planes intersect the input.
// If GenerateFaces is On, it will produce a new closed surface by
// creating new polygonal faces where the input was clipped.
// The GenerateScalars option will add color scalars to the output,
// so that the clipped surfaces can be easily visualized.
// .SECTION Caveats
// This filter does not accept triangle strips as input.
// .SECTION See Also
// vtkOutlineFilter vtkOutlineSource vtkVolumeOutlineSource

#ifndef __vtkClipOutlineWithPlanes_h
#define __vtkClipOutlineWithPlanes_h

#include "vtkPolyDataAlgorithm.h"

class vtkPlaneCollection;
class vtkUnsignedCharArray;
class vtkDoubleArray;
class vtkIdTypeArray;
class vtkCellArray;
class vtkPointData;
class vtkCellData;
class vtkIncrementalPointLocator;
class vtkGenericCell;
class vtkPolygon;
class vtkIdList;

class VTK_EXPORT vtkClipOutlineWithPlanes : public vtkPolyDataAlgorithm
{
public:
  static vtkClipOutlineWithPlanes *New();
  vtkTypeMacro(vtkClipOutlineWithPlanes,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Set the vtkPlaneCollection that holds the clipping planes.
  virtual void SetClippingPlanes(vtkPlaneCollection *planes);
  vtkGetObjectMacro(ClippingPlanes,vtkPlaneCollection);

  // Description:
  // Set whether to add color scalars to the lines.  By default,
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
  // only lines are generated.  Faces will not be generated unless the
  // input has faces, since the output faces are created by clipping
  // the input.
  vtkSetMacro(GenerateFaces, int);
  vtkBooleanMacro(GenerateFaces, int);
  vtkGetMacro(GenerateFaces, int);

  // Description:
  // Set the color for all lines that were part of the original wireframe.
  // The default color is red.  If the the original lines already had
  // color scalars, then this parameter will be ignored.
  vtkSetVector3Macro(BaseColor, double);
  vtkGetVector3Macro(BaseColor, double);

  // Description:
  // Set the color for any new lines that are created as a result of
  // the clipping. The default color is orange.
  vtkSetVector3Macro(ClipColor, double);
  vtkGetVector3Macro(ClipColor, double);

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
  vtkClipOutlineWithPlanes();
  ~vtkClipOutlineWithPlanes();

  vtkPlaneCollection *ClippingPlanes;

  int GenerateScalars;
  int GenerateOutline;
  int GenerateFaces;
  int ActivePlaneId;
  double BaseColor[3];
  double ClipColor[3];
  double ActivePlaneColor[3];

  vtkIncrementalPointLocator *Locator;

  vtkIdList *IdList;
  vtkPolygon *Polygon;

  virtual int ComputePipelineMTime(
    vtkInformation* request, vtkInformationVector** inputVector,
    vtkInformationVector* outputVector, int requestFromOutputPort,
    unsigned long* mtime);

  virtual int RequestData(
    vtkInformation* request, vtkInformationVector** inputVector,
    vtkInformationVector* outputVector);

  void ClipLines(
    vtkPoints *points, vtkDoubleArray *pointScalars,
    vtkIncrementalPointLocator *locator,
    vtkCellArray *inputCells, vtkCellArray *outputLines,
    vtkPointData *inPointData, vtkPointData *outPointData,
    vtkCellData *inCellData, vtkCellData *outLineData);

  void ClipAndContourPolys(
    vtkPoints *points, vtkDoubleArray *pointScalars,
    vtkIncrementalPointLocator *locator,
    vtkCellArray *inputCells, vtkCellArray *outputPolys,
    vtkCellArray *outputLines, vtkPointData *inPointData,
    vtkPointData *outPointData, vtkCellData *inCellData,
    vtkCellData *outPolyData, vtkCellData *outLineData);

  void MakeCutPolys(
    vtkPoints *points, vtkCellArray *lines, vtkIdType firstCell,
    vtkCellArray *polys, double normal[3], vtkCellData *outCD,
    unsigned char colors[3]);

  static void BreakPolylines(
    vtkCellArray *inputlines, vtkCellArray *lines,
    vtkUnsignedCharArray *inputScalars, vtkIdType firstLineScalar,
    vtkUnsignedCharArray *scalars, unsigned char *color);

  static void CopyPolygons(
    vtkCellArray *inputPolys, vtkCellArray *polys,
    vtkUnsignedCharArray *inputScalars, vtkIdType firstPolyScalar,
    vtkUnsignedCharArray *polyScalars, unsigned char colors[3]);

  static void CreateColorValues(
    double color1[3], double color2[3], double color3[3],
    unsigned char colors[3][3]);

private:
  vtkClipOutlineWithPlanes(const vtkClipOutlineWithPlanes&);  // Not implemented.
  void operator=(const vtkClipOutlineWithPlanes&);  // Not implemented.
};

#endif
