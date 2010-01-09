/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkSurfaceCursorShapes.cxx,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkSurfaceCursorShapes.h"
#include "vtkObjectFactory.h"

#include "vtkSurfaceCursor.h"
#include "vtkDataSet.h"
#include "vtkSmartPointer.h"

#include <vtkstd/string>
#include <vtkstd/vector>

vtkCxxRevisionMacro(vtkSurfaceCursorShapes, "$Revision: 1.1 $");
vtkStandardNewMacro(vtkSurfaceCursorShapes);

// A simple container for cursor data
class vtkSurfaceCursorShape
{
public:
  vtkSurfaceCursorShape(const char *name,
                        vtkDataSet *data,
                        int flags) : Name(name), Data(data), Flags(flags) {}; 

  vtkstd::string Name;
  vtkSmartPointer<vtkDataSet> Data;
  int Flags;
};

// A vector of the above, with a VTK-like interface
class vtkSurfaceCursorShapeArray
{
private:
  typedef vtkstd::vector<vtkSurfaceCursorShape> VectorType;
  VectorType Vector;

public:
  static vtkSurfaceCursorShapeArray *New() {
    return new vtkSurfaceCursorShapeArray; };

  void Delete() {
    delete this; };

  void AddItem(const char *name, vtkDataSet *shape, int flags) {
    this->Vector.push_back(vtkSurfaceCursorShape(name, shape, flags)); };

  int GetIndex(const char *name) {
    if (name) {
      for (VectorType::size_type i = 0; i < this->Vector.size(); i++) {
        if (this->Vector[i].Name.compare(name) == 0) {
          return static_cast<int>(i); } } } return -1; };

  const char *GetName(int i) {
    VectorType::size_type j = static_cast<VectorType::size_type>(i);
    if (i < 0 || j >= this->Vector.size()) { return 0; }
    else { return this->Vector[j].Name.c_str(); } };

  vtkDataSet *GetData(int i) {
    VectorType::size_type j = static_cast<VectorType::size_type>(i);
    if (i < 0 || j >= this->Vector.size()) { return 0; }
    else {return this->Vector[j].Data; } };

  int GetFlags(int i) {
    VectorType::size_type j = static_cast<VectorType::size_type>(i);
    if (i < 0 || j >= this->Vector.size()) { return 0; }
    else {return this->Vector[j].Flags; } };
};

//----------------------------------------------------------------------------
vtkSurfaceCursorShapes::vtkSurfaceCursorShapes()
{
  this->NumberOfShapes = 0;
  this->Shapes = vtkSurfaceCursorShapeArray::New();
}

//----------------------------------------------------------------------------
vtkSurfaceCursorShapes::~vtkSurfaceCursorShapes()
{
  this->Shapes->Delete();
}

//----------------------------------------------------------------------------
void vtkSurfaceCursorShapes::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "NumberOfShapes: " << this->NumberOfShapes << "\n";
}

//----------------------------------------------------------------------------
int vtkSurfaceCursorShapes::AddShape(const char *name, vtkDataSet *data,
                                     int flags)
{
  this->Shapes->AddItem(name, data, flags);

  this->NumberOfShapes++;

  return (this->NumberOfShapes - 1);
}

//----------------------------------------------------------------------------
int vtkSurfaceCursorShapes::GetShapeIndex(const char *name)
{
  return this->Shapes->GetIndex(name);
}

//----------------------------------------------------------------------------
const char *vtkSurfaceCursorShapes::GetShapeName(int i)
{
  return this->Shapes->GetName(i);
}

//----------------------------------------------------------------------------
vtkDataSet *vtkSurfaceCursorShapes::GetShapeData(int i)
{
  return this->Shapes->GetData(i);
}

//----------------------------------------------------------------------------
int vtkSurfaceCursorShapes::GetShapeFlags(int i)
{
  return this->Shapes->GetFlags(i);
}

