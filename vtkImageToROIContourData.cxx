/*=========================================================================

  Program:   ToolCursor
  Module:    vtkImageToROIContourData.cxx

  Copyright (c) 2010 David Gobbi
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkImageToROIContourData.h"

#include "vtkROIContourData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkObjectFactory.h"
#include "vtkImageData.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkMath.h"
#include "vtkSynchronizedTemplates2D.h"
#include "vtkStripper.h"

vtkStandardNewMacro(vtkImageToROIContourData);

//----------------------------------------------------------------------------
vtkImageToROIContourData::vtkImageToROIContourData()
{
  this->Value = 0.5;

  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

//----------------------------------------------------------------------------
vtkImageToROIContourData::~vtkImageToROIContourData()
{
}

//----------------------------------------------------------------------------
void vtkImageToROIContourData::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "Value: " << this->Value << "\n";
}

//----------------------------------------------------------------------------
vtkROIContourData* vtkImageToROIContourData::GetOutput()
{
  return vtkROIContourData::SafeDownCast(this->GetOutputDataObject(0));
}

//----------------------------------------------------------------------------
void vtkImageToROIContourData::SetOutput(vtkDataObject* d)
{
  this->GetExecutive()->SetOutputData(0, d);
}

//----------------------------------------------------------------------------
vtkDataObject* vtkImageToROIContourData::GetInput()
{
  return this->GetExecutive()->GetInputData(0, 0);
}

//----------------------------------------------------------------------------
void vtkImageToROIContourData::SetInput(vtkDataObject* input)
{
  vtkAlgorithmOutput *producerPort = 0;

  if (input)
    {
    producerPort = input->GetProducerPort();
    }

  this->SetInputConnection(0, producerPort);
}

//----------------------------------------------------------------------------
int vtkImageToROIContourData::FillOutputPortInformation(
  int, vtkInformation *info)
{
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkROIContourData");
  return 1;
}

//----------------------------------------------------------------------------
int vtkImageToROIContourData::FillInputPortInformation(
  int, vtkInformation *info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
  return 1;
}

//----------------------------------------------------------------------------
int vtkImageToROIContourData::ComputePipelineMTime(
  vtkInformation* vtkNotUsed(request),
  vtkInformationVector** vtkNotUsed(inputVector),
  vtkInformationVector* vtkNotUsed(outputVector),
  int vtkNotUsed(requestFromOutputPort),
  unsigned long* mtime)
{
  unsigned long mTime = this->GetMTime();

  *mtime = mTime;

  return 1;
}

//----------------------------------------------------------------------------
int vtkImageToROIContourData::ProcessRequest(
  vtkInformation* request,
  vtkInformationVector** inputVector,
  vtkInformationVector* outputVector)
{
  // create data object
  if(request->Has(vtkDemandDrivenPipeline::REQUEST_DATA_OBJECT()))
    {
    vtkInformation* info = outputVector->GetInformationObject(0);
    vtkROIContourData *data = vtkROIContourData::SafeDownCast(
      info->Get(vtkDataObject::DATA_OBJECT()));
    if (!data)
      {
      data = vtkROIContourData::New();
      data->SetPipelineInformation(info);
      data->Delete();
      }
    return 1;
    }

  // generate the data
  if(request->Has(vtkDemandDrivenPipeline::REQUEST_DATA()))
    {
    return this->RequestData(request, inputVector, outputVector);
    }

  // tell inputs how to update
  if(request->Has(vtkStreamingDemandDrivenPipeline::REQUEST_UPDATE_EXTENT()))
    {
    vtkInformation* info = inputVector[0]->GetInformationObject(0);
    info->Set(vtkStreamingDemandDrivenPipeline::EXACT_EXTENT(), 1);
    return 1;
    }

  // execute information
  if(request->Has(vtkDemandDrivenPipeline::REQUEST_INFORMATION()))
    {
    return 1;
    }

  return this->Superclass::ProcessRequest(request, inputVector, outputVector);
}

//----------------------------------------------------------------------------
int vtkImageToROIContourData::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // Get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // Get the input and output
  vtkImageData *input = vtkImageData::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkROIContourData *output = vtkROIContourData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  // Create filters to do the work
  vtkSynchronizedTemplates2D *synctemp =
    vtkSynchronizedTemplates2D::New();
  vtkStripper *stripper = vtkStripper::New();
  stripper->SetInputConnection(synctemp->GetOutputPort());

  // Create an image to use as input to SyncTemp
  vtkImageData *slice = vtkImageData::New();
  slice->SetOrigin(input->GetOrigin());
  slice->SetSpacing(input->GetSpacing());
  slice->SetScalarType(input->GetScalarType());
  slice->SetNumberOfScalarComponents(input->GetNumberOfScalarComponents());
  slice->AllocateScalars();
  synctemp->SetInput(slice);
  synctemp->SetValue(0, this->Value);

  // Go through the input slice by slice
  int extent[6];
  input->GetExtent(extent);
  int zMin = extent[4];
  int zMax = extent[5];
  vtkIdType sliceSize = extent[1] - extent[0] + 1;
  sliceSize *= extent[3] - extent[2] + 1;
  sliceSize *= input->GetNumberOfScalarComponents();

  for (int zIdx = zMin; zIdx <= zMax; zIdx++)
    {
    // Set extent to just one slice
    extent[4] = zIdx;
    extent[5] = zIdx;
    slice->SetWholeExtent(extent);
    slice->SetExtent(extent);
    slice->GetPointData()->GetScalars()->SetVoidArray(
      input->GetScalarPointerForExtent(extent), sliceSize, 1);

    // Process and get output
    stripper->Update();
    vtkPolyData *sliceContours = stripper->GetOutput();
    vtkPoints *slicePoints = sliceContours->GetPoints();
    vtkCellArray *sliceLines = sliceContours->GetLines();

    // Add contours to output
    vtkIdType numCells = sliceLines->GetNumberOfCells();
    int contourId = output->GetNumberOfContours();
    output->SetNumberOfContours(contourId + numCells);
    vtkIdType loc = 0;
    for (vtkIdType j = 0; j < numCells; j++)
      {
      vtkIdType numPts, *ptIds;
      sliceLines->GetCell(loc, numPts, ptIds);
      loc += numPts + 1;
      vtkPoints *points = vtkPoints::New();
      vtkIdType n = numPts;
      bool closed = (ptIds[0] == ptIds[numPts-1]);
      n -= closed;
      points->SetNumberOfPoints(n);
      for (int i = 0; i < n; i++)
        {
        double p[3];
        slicePoints->GetPoint(ptIds[i], p);
        points->SetPoint(i, p);
        }
      output->SetContourPoints(contourId, points);
      output->SetContourType(contourId,
         (closed ? vtkROIContourData::CLOSED_PLANAR :
                   vtkROIContourData::OPEN_PLANAR));
      points->Delete();
      contourId++;
      }
    }

  // Free temporary objects
  synctemp->SetInput(0);
  stripper->SetInput(0);
  synctemp->Delete();
  stripper->Delete();
  slice->Delete();

  return 1;
}
