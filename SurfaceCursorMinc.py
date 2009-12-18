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

#---------------------------------------------------------
# the cone points along the -x axis
coneSource = vtk.vtkConeSource()
coneSource.CappingOn()
coneSource.SetHeight(12)
coneSource.SetRadius(5)
coneSource.SetResolution(31)
coneSource.SetCenter(6,0,0)
coneSource.SetDirection(-1,0,0)

coneMapper = vtk.vtkDataSetMapper()
coneMapper.SetInputConnection(coneSource.GetOutputPort())

redCone = vtk.vtkActor()
redCone.PickableOff()
redCone.SetMapper(coneMapper)
redCone.GetProperty().SetColor(1,0,0)

greenCone = vtk.vtkActor()
greenCone.PickableOff()
greenCone.SetMapper(coneMapper)
greenCone.GetProperty().SetColor(0,1,0)

#ren.AddViewProp(redCone)
#ren.AddViewProp(greenCone)

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

# A function to point an actor along a vector
def PointCone(actor,nx,ny,nz): #,px,py,pz):
    actor.SetOrientation(0.0, 0.0, 0.0)
    n = math.sqrt(nx**2 + ny**2 + nz**2)
    if (nx < 0.0):
        actor.RotateWXYZ(180, 0, 1, 0)
        actor.RotateWXYZ(180, (nx-n)*0.5, ny*0.5, nz*0.5)
    else:
        actor.RotateWXYZ(180, (nx+n)*0.5, ny*0.5, nz*0.5)

def MoveCursor(iren,event=""):
    level = 0
    if iren.GetShiftKey():
      level = level | 1
    if iren.GetControlKey():
      level = level | 2
    if event == "KeyPressEvent":
      if iren.GetKeySym() == "Control_L":
        level = level | 2
      elif iren.GetKeySym() == "Shift_L":
        level = level | 1
      print "KeyPress shift=%i, control=%i, keysym=\"%s\"" % (iren.GetShiftKey(), iren.GetControlKey(), iren.GetKeySym())
    elif event == "KeyReleaseEvent":
      if iren.GetKeySym() == "Control_L":
        level = level & ~2
      elif iren.GetKeySym() == "Shift_L":
        level = level & ~1
      print "KeyRelease shift=%i, control=%i, keysym=\"%s\"" % (iren.GetShiftKey(), iren.GetControlKey(), iren.GetKeySym())
    cursor.SetLevel(level)
    iren.Render()

def EnterRenWin(iren, event=""):
    renWin.HideCursor()

def LeaveRenWin(iren, event=""):
    renWin.ShowCursor()

def OnRender(ren,event=""):
    x,y = iren.GetEventPosition()
    cursor.SetDisplayPosition(x,y)
    picker = cursor.GetPicker()
    p = picker.GetPickPosition()
    n = picker.GetPickNormal()
    redCone.SetPosition(p[0],p[1],p[2])
    PointCone(redCone,n[0],n[1],n[2])
    greenCone.SetPosition(p[0],p[1],p[2])
    PointCone(greenCone,-n[0],-n[1],-n[2])

#---------------------------------------------------------
# custom interaction
iren.AddObserver("MouseMoveEvent", MoveCursor)
iren.AddObserver("KeyPressEvent", MoveCursor)
iren.AddObserver("KeyReleaseEvent", MoveCursor)
iren.AddObserver("EnterEvent", EnterRenWin)
iren.AddObserver("LeaveEvent", LeaveRenWin)
ren.AddObserver("StartEvent", OnRender)

iren.Start()
