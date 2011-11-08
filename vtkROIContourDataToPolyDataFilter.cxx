/*=========================================================================

  Program:   ToolCursor
  Module:    vtkROIContourDataToPolyDataFilter.cxx

  Copyright (c) 2010 David Gobbi
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkROIContourDataToPolyDataFilter.h"

#include "vtkROIContourData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkImageData.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkPlane.h"
#include "vtkMath.h"

vtkStandardNewMacro(vtkROIContourDataToPolyDataFilter);

vtkCxxSetObjectMacro(vtkROIContourDataToPolyDataFilter,SelectionPlane,vtkPlane);

//----------------------------------------------------------------------------
vtkROIContourDataToPolyDataFilter::vtkROIContourDataToPolyDataFilter()
{
  this->SelectionPlane = NULL;
  this->SelectionPlaneTolerance = 0.5;
}

//----------------------------------------------------------------------------
vtkROIContourDataToPolyDataFilter::~vtkROIContourDataToPolyDataFilter()
{
  if (this->SelectionPlane)
    {
    this->SelectionPlane->Delete();
    }
}

//----------------------------------------------------------------------------
void vtkROIContourDataToPolyDataFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "SelectionPlane: ";
  if (this->SelectionPlane)
    {
    os << this->SelectionPlane << "\n";
    }
  else
    {
    os << "(none)\n";
    }
  os << indent << "SelectionPlaneTolerance: "
     << this->SelectionPlaneTolerance << "\n";
}

//----------------------------------------------------------------------------
int vtkROIContourDataToPolyDataFilter::FillInputPortInformation(
  int, vtkInformation *info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkROIContourData");
  return 1;
}

//----------------------------------------------------------------------------
int vtkROIContourDataToPolyDataFilter::ComputePipelineMTime(
  vtkInformation* vtkNotUsed(request),
  vtkInformationVector** vtkNotUsed(inputVector),
  vtkInformationVector* vtkNotUsed(outputVector),
  int vtkNotUsed(requestFromOutputPort),
  unsigned long* mtime)
{
  unsigned long mTime = this->GetMTime();

  vtkPlane *plane = this->SelectionPlane;

  if (plane)
    {
    unsigned long planeMTime = plane->GetMTime();
    if (planeMTime > mTime)
      {
      mTime = planeMTime;
      }
    }

  *mtime = mTime;

  return 1;
}

//----------------------------------------------------------------------------
int vtkROIContourDataToPolyDataFilter::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // Get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // Get the input and output
  vtkROIContourData *input = vtkROIContourData::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData *output = vtkPolyData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  // For selecting which contours to include in the output
  vtkPlane *plane = this->SelectionPlane;
  double tol = this->SelectionPlaneTolerance;

  // The output points
  vtkPoints *outPoints = vtkPoints::New(VTK_DOUBLE);
  vtkCellArray *lines = 0;
  vtkCellArray *verts = 0;

  // Go through all the contours
  int n = input->GetNumberOfContours();
  for (int i = 0; i < n; i++)
    {
    vtkPoints *points = input->GetContourPoints(i);
    int t = input->GetContourType(i);

    if (points)
      {
      vtkIdType m = points->GetNumberOfPoints();
      bool includeContour = true;

      if (plane)
        {
        for (int j = 0; j < m; j++)
          {
          double p[3];
          points->GetPoint(j, p);
          double d = plane->DistanceToPlane(p);
          if (d < -tol || d > tol)
            {
            includeContour = false;
            break;
            }
          }
        }

      // Include this contour in the output
      if (includeContour)
        {
        vtkCellArray *cells = 0;
        vtkIdType cellSize = m;
        if (t == vtkROIContourData::POINT)
          {
          if (!verts)
            {
            verts = vtkCellArray::New();
            cells = verts;
            }
          }
        else
          {
          if (!lines)
            {
            lines = vtkCellArray::New();
            cells = lines;
            }
          if (t == vtkROIContourData::CLOSED_PLANAR)
            {
            cellSize = m + 1;
            }
          }

        // Generate the cell
        cells->InsertNextCell(cellSize);
        vtkIdType firstPointId = outPoints->GetNumberOfPoints();

        for (int j = 0; j < m; j++)
          {
          double p[3];
          points->GetPoint(j, p);
          vtkIdType pointId = outPoints->InsertNextPoint(p);
          cells->InsertCellPoint(pointId);
          }

        // Close the contour, if necessary
        if (cellSize > m)
          {
          cells->InsertCellPoint(firstPointId);
          }
        }
      }
    }

  output->SetPoints(outPoints);
  output->SetLines(lines);
  output->SetVerts(verts);

  if (outPoints)
    {
    outPoints->Delete();
    }
  if (lines)
    {
    lines->Delete();
    }
  if (verts)
    {
    verts->Delete();
    }

  return 1;
}
