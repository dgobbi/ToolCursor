#!/usr/bin/env python

import math
import sys
import vtk
import surfacecursor

from vtk.util.misc import vtkGetDataRoot
VTK_DATA_ROOT = vtkGetDataRoot()

if len(sys.argv) > 1:
  filename = sys.argv[1]
else:
  filename = "/beck/data/dgobbi/Atamai/atamai/examples/sbrain.mnc"

#---------------------------------------------------------
# renderer and interactor
ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren)
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)

#---------------------------------------------------------
# read the volume
reader1 = vtk.vtkMINCImageReader()
reader1.SetFileName(filename)
#reader = vtk.vtkImageReader2()
#reader.SetDataExtent(0,63,0,63,0,92)
#reader.SetFileNameSliceOffset(1)
#reader.SetDataScalarTypeToUnsignedShort()
#reader.SetDataByteOrderToLittleEndian()
#reader.SetFilePrefix(str(VTK_DATA_ROOT) + "/Data/headsq/quarter")
#reader.SetDataSpacing(3.2,3.2,1.5)

reader = vtk.vtkImageShiftScale()
reader.SetInput(reader1.GetOutput())
reader.SetOutputScalarTypeToUnsignedShort()

#---------------------------------------------------------
# set up the volume rendering

volumeMapper = vtk.vtkVolumeTextureMapper3D()
volumeMapper.SetInput(reader.GetOutput())
#volumeFunction = vtk.vtkVolumeRayCastCompositeFunction()
#volumeMapper.SetVolumeRayCastFunction(volumeFunction)
volumeMapper.CroppingOn()
# (-90.0, 90.0, -126.0, 90.0, -72.0, 108.0)
volumeMapper.SetCroppingRegionPlanes((0.0, 90.0, -126.0, 0.0, -72.0, 0.0))
volumeMapper.SetCroppingRegionFlagsToFence()

volumeColor = vtk.vtkColorTransferFunction()
volumeColor.AddRGBPoint(0,0.0,0.0,0.0)
volumeColor.AddRGBPoint(180,0.3,0.1,0.2)
volumeColor.AddRGBPoint(1200,1.0,0.7,0.6)
volumeColor.AddRGBPoint(2500,1.0,1.0,0.9)

volumeScalarOpacity = vtk.vtkPiecewiseFunction()
volumeScalarOpacity.AddPoint(0,0.0)
volumeScalarOpacity.AddPoint(180,0.0)
volumeScalarOpacity.AddPoint(1200,0.2)
volumeScalarOpacity.AddPoint(2500,0.8)

volumeGradientOpacity = vtk.vtkPiecewiseFunction()
volumeGradientOpacity.AddPoint(0,0.0)
volumeGradientOpacity.AddPoint(90,0.5)
volumeGradientOpacity.AddPoint(100,1.0)

volumeProperty = vtk.vtkVolumeProperty()
volumeProperty.SetColor(volumeColor)
volumeProperty.SetScalarOpacity(volumeScalarOpacity)
#volumeProperty.SetGradientOpacity(volumeGradientOpacity)
volumeProperty.SetInterpolationTypeToLinear()
volumeProperty.ShadeOff()
volumeProperty.SetAmbient(0.6)
volumeProperty.SetDiffuse(0.6)
volumeProperty.SetSpecular(0.1)

volume = vtk.vtkVolume()
volume.SetMapper(volumeMapper)
volume.SetProperty(volumeProperty)

#---------------------------------------------------------
# Do the surface rendering
skinExtractor = vtk.vtkMarchingCubes()
skinExtractor.SetInputConnection(reader.GetOutputPort())
skinExtractor.SetValue(0,500)

skinNormals = vtk.vtkPolyDataNormals()
skinNormals.SetInputConnection(skinExtractor.GetOutputPort())
skinNormals.SetFeatureAngle(60.0)

skinStripper = vtk.vtkStripper()
skinStripper.SetMaximumLength(10)
skinStripper.SetInputConnection(skinNormals.GetOutputPort())

skinLocator = vtk.vtkCellLocator()
skinLocator.SetDataSet(skinStripper.GetOutput())
skinLocator.LazyEvaluationOn()

skinMapper = vtk.vtkPolyDataMapper()
skinMapper.SetInputConnection(skinStripper.GetOutputPort())
skinMapper.ScalarVisibilityOff()

skinProperty = vtk.vtkProperty()
skinProperty.SetColor(1.0,1.0,0.9)

skin = vtk.vtkActor()
skin.SetMapper(skinMapper)
skin.SetProperty(skinProperty)

#---------------------------------------------------------
# Create an image actor
table = vtk.vtkLookupTable()
table.SetRange(0,2000)
table.SetRampToLinear()
table.SetValueRange(0,1)
table.SetHueRange(0,0)
table.SetSaturationRange(0,0)

mapToColors = vtk.vtkImageMapToColors()
mapToColors.SetInputConnection(reader.GetOutputPort())
mapToColors.SetLookupTable(table)
mapToColors.GetOutput().Update()

imageActor = vtk.vtkImageActor()
imageActor.SetInput(mapToColors.GetOutput())
imageActor.SetDisplayExtent(32,32,0,63,0,92)

#---------------------------------------------------------
# make a transform and some clipping planes
transform = vtk.vtkTransform()
transform.RotateWXYZ(-20,0.0,-0.7,0.7)

#volume.SetUserTransform(transform)
#skin.SetUserTransform(transform)
#imageActor.SetUserTransform(transform)

c = volume.GetCenter()

volumeClip = vtk.vtkPlane()
volumeClip.SetNormal(0,1,0)
volumeClip.SetOrigin(c[0],c[1],c[2])

skinClip = vtk.vtkPlane()
skinClip.SetNormal(0,0,1)
skinClip.SetOrigin(c[0],c[1],c[2])

#volumeMapper.AddClippingPlane(volumeClip)
skinMapper.AddClippingPlane(skinClip)

#---------------------------------------------------------
ren.AddViewProp(volume)
#ren.AddViewProp(skin)
#ren.AddViewProp(imageActor)

camera = ren.GetActiveCamera()
camera.SetFocalPoint(c[0],c[1],c[2])
camera.SetPosition(c[0] - 500,c[1] + 100,c[2] + 100)
camera.SetViewUp(0,0,1)

renWin.Render()

#ren.SetBackground(0.5, 0.5, 0.5)
ren.ResetCameraClippingRange()
#---------------------------------------------------------
# the picker
picker = vtk.vtkSurfacePicker()
picker.SetTolerance(1e-6)
picker.SetVolumeOpacityIsovalue(0.1)
picker.AddLocator(skinLocator)

cursor = vtk.vtkSurfaceCursor()
cursor.SetRenderer(ren)
cursor.SetScale(1)
cursor.BindInteractor(iren)

def OnRender(ren,event=""):
    x,y = iren.GetEventPosition()
    cursor.SetDisplayPosition(x,y)

#---------------------------------------------------------
# custom interaction
ren.AddObserver("StartEvent", OnRender)

iren.Start()
