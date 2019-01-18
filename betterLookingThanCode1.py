# -*- coding: utf-8 -*-
"""
Created on Thu Jan 17 15:24:04 2019

@author: mhmdk
"""

import opensim as osim
import matplotlib.pyplot as plt
import re
from tkinter import*

#load the model we want to use from corresponding directory
newModel = osim.Model("C:/OpenSim 3.3/Models/Tug_of_War/Tug_of_War_Millard.osim")
newMuscleBase = newModel.getMuscles().get(0)
newMuscle = osim.Millard2012EquilibriumMuscle.safeDownCast(newMuscleBase)
initState = newModel.initSystem()

alpha = DoubleVar()
act = DoubleVar()
deAct = DoubleVar()
ut = DoubleVar()
Ao = DoubleVar()
vm = 1e-4
fo = DoubleVar()
lo = DoubleVar()
ls = DoubleVar()
lMT = float(ls.get())+ float(lo.get())
vMax = DoubleVar()

def alphaCom(num):
    newMuscle.set_pennation_angle_at_optimal(float(num))
    #print("angle: " + str(newMuscle.get_pennation_angle_at_optimal()))

def actCom(num):
    diff = deAct.get() - float(num)
    if diff < 0.03:
        setter = deAct.get() + (0.03 - diff)
        scaleDeAct.set(setter)
    if diff > 0.04:
        setter = deAct.get() - (diff - 0.04)
        scaleDeAct.set(setter)
    newMuscle.set_activation_time_constant(float(num))
    #print("activation time constant: " + str(newMuscle.get_activation_time_constant()))
    
def deActCom(num):
    diff = float(num) - act.get()
    if diff < 0.03:
        setter = act.get() - (0.03 - diff)
        scaleAct.set(setter)
    if diff > 0.04:
        setter = act.get() + (diff - 0.04)
        scaleAct.set(setter)
    newMuscle.set_deactivation_time_constant(float(num))
    #print("deactivation time constant: " + str(newMuscle.get_deactivation_time_constant()))
    
def areaCom(num):
    area = float(num)
    fibLen = lo.get()
    vmtest = area*fibLen
    if vmtest > 1e-4:
        fibTester = 1e-4/area
        if (fibTester >= 0.05) & (fibTester <= 0.2):
            scaleLo.set(fibTester)
        else:
            if fibTester > 0.2:
                scaleLo.set(0.2)
            else:
                scaleLo.set(0.05)
    force = 350000*area
    scaleFo.set(force)
    # idk what to put for the transferMuscleProperties for area

def forceCom(num):
    force = float(num)
    fibLen = lo.get()
    maxV = vMax.get()
    power = force*fibLen*maxV
    if power > 175.0:
        vTester = 175.0/(force*fibLen)
        if (vTester >= 2.0) & (vTester <= 10.0):
            scaleVmax.set(vTester)
        else:
            scaleVmax.set(2)
            fibLenSetter = 175/(force*2)
            scaleLo.set(fibLenSetter)
    area = force/350000
    scaleAo.set(area)
    newMuscle.set_max_isometric_force(float(num))
    #print("maximum isometric force: " + str(newMuscle.get_max_isometric_force()))

def maxVelCom(num):
    vel = float(num)
    force = fo.get()
    fibLen = lo.get()
    if vel*fibLen*force > 175:
        fTester = 175/(fibLen*vel)
        if (fTester < 1750) & (fTester > 87.5):
            scaleFo.set(fTester)
        else:
            force = 87.5
            scaleFo.set(87.5)
            scaleLo.set(175/(force*vel))
    newMuscle.set_max_contraction_velocity(float(num))
    #print("maximum contraction velocity: " + str(newMuscle.get_max_contraction_velocity()))

def fibLenCom(num):
    fibLen = float(num)
    slackLen = ls.get()
    force = fo.get()
    maxV = vMax.get()
    area = Ao.get()

    if force*maxV*fibLen > 175:
        fTester = 175/(fibLen*maxV)
        if (fTester > 87.5) & (fTester < 1750):
            scaleFo.set(fTester)
        else:
            scaleFo.set(87.5)
            maxV = 175/(87.5*fibLen)
            scaleVmax.set(maxV)
    sum  =  fibLen + slackLen
    if sum > 0.45:
        diff = sum - 0.45
        setter = slackLen - diff
        scaleLs.set(setter)
    elif sum < 0.15:
        diff  =  0.15 - sum
        setter = slackLen + diff
        scaleLs.set(setter)

    if area*fibLen > 1e-4:
        setter = 1e-4/fibLen
        scaleAo.set(setter)
    setter = float(lo.get()) + float(ls.get())
    labelTotLen.configure(text="distance between origin and attachment: %0.2f" % setter)
    newMuscle.set_optimal_fiber_length(float(num))
    blockLoc = newMuscle.getGeometryPath().getPathPointSet().get(1).getLocation().get(2)
    # Left muscle
    if blockLoc > 0:
		groundLoc = blockLoc + fibLen + slackLen
	# Right muscle
    else:
        groundLoc = blockLoc - fibLen - slackLen
	vector = osim.Vec3(0, 0.5, groundLoc)
	newMuscle.updGeometryPath().getPathPointSet().get(0).setLocation(initState,vector)
    #print("optimal fiber length: " + str(newMuscle.get_optimal_fiber_length()))

def slackLenCom(num):
    slackLen = float(num)
    fibLen = lo.get()
    totLen = fibLen + slackLen
    if totLen > 0.45:
        diff = totLen - 0.45
        setter = fibLen - diff
        scaleLo.set(setter)
    elif totLen < 0.15:
        diff = 0.15 - totLen
        setter = diff + fibLen
        scaleLo.set(setter)
    setter = float(lo.get()) + float(ls.get())
    labelTotLen.configure(text = "distance between origin and attachment: %0.2f" %setter)
    newMuscle.set_tendon_slack_length(float(num))
    blockLoc = newMuscle.getGeometryPath().getPathPointSet().get(1).getLocation().get(2)
    # Left muscle
    if blockLoc > 0:
		groundLoc = blockLoc + fibLen + slackLen
	# Right muscle
    else:
        groundLoc = blockLoc - fibLen - slackLen
	vector = osim.Vec3(0, 0.5, groundLoc)
	newMuscle.updGeometryPath().getPathPointSet().get(0).setLocation(initState, vector)
    #print("tendon slack length: " + str(newMuscle.get_tendon_slack_length()))

def main():
    Label(root, text = "Angle in degrees: ").pack()
    scaleAlpha = Scale(root, variable = alpha, orient = HORIZONTAL, from_ = 0, to = 30, command = alphaCom)
    scaleAlpha.pack()
    scaleAlpha.set(15)
    
    Label(root, text = "Total excitation: ").pack()
    scaleUt = Scale(root, variable = ut, orient = HORIZONTAL,  from_ = 0, to = 0.5, resolution = 0.05)
    scaleUt.pack()
    scaleUt = 0.25
    
    actT = Label(root, text = "Activation time constant: ")
    actT.pack()
    scaleAct = Scale(root, variable = act, orient = HORIZONTAL,  from_ = 0.01, to = 0.02, resolution = 0.0005, command = actCom)
    scaleAct.pack()
    scaleAct.set(0.015)
    
    deActT = Label(root, text = "Deactivation time constant: ")
    deActT.pack()
    scaleDeAct = Scale(root, variable = deAct, orient = HORIZONTAL, from_ = 0.04, to = 0.06, resolution = 0.0005, command = deActCom)
    scaleDeAct.pack()
    scaleDeAct.set(0.05)
    
    #the constraint does not describe an upper limit for the area
    AoT = Label(root, text = "Area of muscle: ").pack()
    scaleAo = Scale(root, variable = Ao, orient = HORIZONTAL, from_ = 0.00025, to = 0.005, resolution = 0.00005, command = areaCom)
    scaleAo.pack()
    scaleAo.set(0.002625)
    
    #also lower and upper limits are not defined
    forceT = Label(root, text = "Maximum isometric muscle fiber force: ").pack()
    scaleFo = Scale(root, variable = fo, orient = HORIZONTAL, from_ = 87.5, to = 1750, resolution = 0.5, command = forceCom)
    scaleFo.pack()
    scaleFo.set(918)
    
    loT = Label(root, text = "Optimal muscle fiber: ")
    loT.pack()
    scaleLo = Scale(root, variable = lo, orient = HORIZONTAL, from_ = 0.05, to = 0.2, resolution = 0.01, command = fibLenCom)
    scaleLo.pack()
    scaleLo.set(0.1)
    
    # upper limit is not specified
    slackT = Label(root, text = "Tendon slack length: ").pack()
    scaleLs = Scale(root, variable = ls, orient = HORIZONTAL, from_ = 0.1, to = 0.4, resolution = 0.01, command = slackLenCom)
    scaleLs.pack()
    scaleLs.set(0.25)
    
    vMaxT = Label(root, text = "Maximum contraction velocity: ").pack()
    scaleVmax = Scale(root, variable = vMax, orient = HORIZONTAL, from_ = 2, to = 10, command = maxVelCom)
    scaleVmax.pack()
    scaleVmax.set(6)
    
    labelVm = Label(root, text = "Muscle Volume (cm^3): " + str(vm*1e6))
    labelVm.pack()
        