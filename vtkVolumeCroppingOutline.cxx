/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkVolumeCroppingOutline.cxx,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkVolumeCroppingOutline.h"

#include "vtkDataSet.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkImageData.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkUnsignedCharArray.h"
#include "vtkVolumeMapper.h"
#include "vtkMath.h"

vtkCxxRevisionMacro(vtkVolumeCroppingOutline, "$Revision: 1.2 $");
vtkStandardNewMacro(vtkVolumeCroppingOutline);

vtkCxxSetObjectMacro(vtkVolumeCroppingOutline,VolumeMapper,vtkVolumeMapper);

//----------------------------------------------------------------------------
vtkVolumeCroppingOutline::vtkVolumeCroppingOutline ()
{
  this->VolumeMapper = 0;
  this->UseColorScalars = 0;
  this->ActivePlane = -1;

  this->Color[0] = 1.0;
  this->Color[1] = 0.0;
  this->Color[2] = 0.0;

  this->ActivePlaneColor[0] = 1.0;
  this->ActivePlaneColor[1] = 1.0;
  this->ActivePlaneColor[2] = 0.0;

  this->SetNumberOfInputPorts(0);
}

//----------------------------------------------------------------------------
vtkVolumeCroppingOutline::~vtkVolumeCroppingOutline ()
{
  if (this->VolumeMapper)
    {
    this->VolumeMapper->Delete();
    this->VolumeMapper = 0;
    }
}

//----------------------------------------------------------------------------
void vtkVolumeCroppingOutline::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "VolumeMapper: ";
  if (this->VolumeMapper)
    {
    os << this->VolumeMapper << "\n";
    }
  else
    {
    os << "(none)\n";
    } 

  os << indent << "UseColorScalars: "
     << (this->UseColorScalars ? "On\n" : "Off\n" );

  os << indent << "Color: " << this->Color[0] << ", "
     << this->Color[1] << ", " << this->Color[2] << "\n";

  os << indent << "ActivePlane: " << this->ActivePlane << "\n";

  os << indent << "ActivePlaneColor: " << this->ActivePlaneColor[0] << ", "
     << this->ActivePlaneColor[1] << ", " << this->ActivePlaneColor[2] << "\n";
}

//----------------------------------------------------------------------------
int vtkVolumeCroppingOutline::ComputeCubePlanes(double planes[3][4])
{
  for (int i = 0; i < 3; i++)
    {
    int j0 = 2*i;
    int j1 = 2*i + 1;

    double a = this->Bounds[j0];
    double b = this->CroppingRegionPlanes[j0];
    double c = this->CroppingRegionPlanes[j1];
    double d = this->Bounds[j1];

    if (a > d || b > c)
      {
      return 0;
      }

    if (b < a) { b = a; };
    if (b > d) { b = d; };
    if (c < a) { c = a; };
    if (c > d) { c = d; };

    planes[i][0] = a;
    planes[i][1] = b;
    planes[i][2] = c;
    planes[i][3] = d;
    }
}

//----------------------------------------------------------------------------
int vtkVolumeCroppingOutline::ComputePipelineMTime(
  vtkInformation* vtkNotUsed(request),
  vtkInformationVector** vtkNotUsed(inputVector),
  vtkInformationVector* vtkNotUsed(outputVector),
  int vtkNotUsed(requestFromOutputPort),
  unsigned long* mtime)
{
  unsigned long mTime = this->GetMTime();
  if (this->VolumeMapper)
    {
    int mapperMTime = this->VolumeMapper->GetMTime();
    if (mapperMTime > mTime)
      {
      mTime = mapperMTime;
      }
    vtkImageData *input = this->VolumeMapper->GetInput();
    if (input)
      {
      // Need to do this because we are not formally connected
      // to the Mapper's pipeline
      input->UpdateInformation();
      unsigned long pipelineMTime = input->GetPipelineMTime();
      if (pipelineMTime > mTime)
        {
        mTime = pipelineMTime;
        }
      }
    }

  *mtime = mTime;

  return 1;
}

//----------------------------------------------------------------------------
int vtkVolumeCroppingOutline::RequestInformation(
  vtkInformation* vtkNotUsed(request),
  vtkInformationVector** vtkNotUsed(inputVector),
  vtkInformationVector* vtkNotUsed(outputVector))
{
  // Get the mapper's input, since this is the most convenient
  // place to do so.

  if (!this->VolumeMapper)
    {
    vtkErrorMacro("No VolumeMapper has been set.");
    return 0;
    }

  this->Cropping = this->VolumeMapper->GetCropping();
  this->CroppingRegionFlags = this->VolumeMapper->GetCroppingRegionFlags();
  this->VolumeMapper->GetCroppingRegionPlanes(this->CroppingRegionPlanes);

  vtkImageData *data = this->VolumeMapper->GetInput();

  if (!data)
    {
    vtkErrorMacro("The VolumeMapper does not have an input set.");
    return 0;
    }

  // Don't have to update mapper's input, since it was done in
  // ComputePipelineMTime.
  // data->UpdateInformation();

  // Don't call GetBounds because we need whole extent

  double spacing[3];
  double origin[3];
  int extent[6];

  data->GetSpacing(spacing);
  data->GetOrigin(origin);
  data->GetWholeExtent(extent);

  for (int i = 0; i < 3; i++)
    {
    int j0 = 2*i;
    int j1 = j0+1;

    if (extent[j0] > extent[j1])
      {
      vtkMath::UninitializeBounds(this->Bounds);
      break;
      }

    if (spacing[i] > 0)
      {
      this->Bounds[j0] = origin[i] + spacing[i]*extent[j0];
      this->Bounds[j1] = origin[i] + spacing[i]*extent[j1];
      }
    else
      {
      this->Bounds[j0] = origin[i] + spacing[i]*extent[j1];
      this->Bounds[j1] = origin[i] + spacing[i]*extent[j0];
      }
    }

  return 1;
}

//----------------------------------------------------------------------------
int vtkVolumeCroppingOutline::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get the info object
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the output
  vtkPolyData *output = vtkPolyData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkDebugMacro(<< "Creating cropping region outline");

  double planes[3][4];

  if (!this->ComputeCubePlanes(planes))
    {
    }

  // Add an offset to the planes here

  int flags = this->CroppingRegionFlags;

  // The points
  vtkPoints *points = vtkPoints::New();
  points->Allocate(64);

  // Create the array of 64 points
  int ptId = 0;
  for (int i = 0; i < 4; i++)
    {
    for (int j = 0; j < 4; j++)
      {
      for (int k = 0; k < 4; k++)
        {
        points->InsertPoint(ptId++, planes[0][k], planes[1][j], planes[2][i]);
        }
      }
    }

  // The active plane, which gets the second color
  int activePlane = this->ActivePlane;
  if (activePlane > 5) { activePlane = -1; };

  // The colors
  unsigned char colors[2][3];
  colors[0][0] = colors[0][1] = colors[0][2] = 255;
  colors[1][0] = colors[1][1] = colors[1][2] = 255;

  // The scalars to color the lines
  vtkUnsignedCharArray *scalars = 0;
  if (this->UseColorScalars)
    {
    // Convert the two colors to unsigned char
    double *dcolors[2];
    dcolors[0] = this->Color;
    dcolors[1] = this->ActivePlaneColor;

    for (int i = 0; i < 2; i++)
      {
      for (int j = 0; j < 3; j++)
        {
        double val = dcolors[i][j];
        if (val < 0) { val = 0; }
        if (val > 1) { val = 1; }
        colors[i][j] = static_cast<unsigned char>(val*255);
        }
      }

    scalars = vtkUnsignedCharArray::New();
    scalars->SetNumberOfComponents(3);
    }

  // The lines
  vtkCellArray *lines = vtkCellArray::New();

  // Loop over the three dimensions and create the lines
  for (int dim0 = 0; dim0 < 3; dim0++)
    {
    // Compute the other two dimension indices
    int dim1 = (dim0+1)%3;
    int dim2 = (dim0+2)%3;

    // Indices into the cubes
    int idx[3];

    // Loop over dimensions other than "dim"
    for (int i = 0; i < 4; i++)
      {
      idx[dim2] = i;

      for (int j = 0; j < 4; j++)
        {
        idx[dim1] = j;

        // Loop over line segments
        for (int k = 0; k < 3; k++)
          {
          idx[dim0] = k;

          // The first point in the segment, and the increment to the next
          int pointId = idx[2]*16 + idx[1]*4 + idx[0];
          int pointInc = (1 << (2*dim0));

          // Loop through the four cubes adjacent to the line segment
          int bitCheck = 0;
          int cidx[3];
          cidx[dim0] = idx[dim0];
          for (int ii = 0; ii < 2; ii++)
            {
            // First get idx[dim1]-1, then idx[dim1]
            cidx[dim1] = idx[dim1] + ii - 1;
            for (int jj = 0; jj < 2; jj++)
              {
              // First get idx[dim2]-1, then idx[dim2], but reverse
              // the order when ii loop is on its second iteration
              cidx[dim2] = idx[dim2] + (ii^jj) - 1;
              int flagval = 0;
              if (cidx[dim1] >= 0 && cidx[dim1] < 3 &&
                  cidx[dim2] >= 0 && cidx[dim2] < 3)
                {
                int flagbit = cidx[2]*9 + cidx[1]*3 + cidx[0];
                flagval = ((flags >> flagbit) & 1);
                }
              bitCheck <<= 1;
              bitCheck |= flagval;
              }
            }

          // Whether we need a line depends on the the value of bitCheck.
          // Values 0000, 0011, 0110, 1100, 1001, 1111 don't need lines. 
          // Build a bitfield rather than a lookup table, it's faster.
          const int noLineValues = ((1 << 0x0) | (1 << 0x3) | (1 << 0x6) |
                                    (1 << 0x9) | (1 << 0xc) | (1 << 0xf)); 

          if (((noLineValues >> bitCheck) & 1) == 0)
            {
            lines->InsertNextCell(2);
            lines->InsertCellPoint(pointId);
            lines->InsertCellPoint(pointId + pointInc);
            if (scalars)
              {
              unsigned char *color = colors[0];
              if (activePlane >= 0)
                {
                int planeDim = (activePlane >> 1); // same as "/ 2"
                int planeIdx = 1 + (activePlane & 1); // same as "% 2"
                if ((planeDim == dim2 && i == planeIdx) ||
                    (planeDim == dim1 && j == planeIdx))
                  {
                  color = colors[1];
                  }
                }  
              scalars->InsertNextTupleValue(color);
              }
            }
          }
        }
      }
    }

  output->SetPoints(points);
  points->Delete();

  output->SetLines(lines);
  lines->Delete();

  output->GetCellData()->SetScalars(scalars);
  if (scalars)
    {
    scalars->Delete();
    }

  return 1;
}

