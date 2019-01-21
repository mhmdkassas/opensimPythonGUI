# ----------------------------------------------------------------------- #
# The OpenSim API is a toolkit for musculoskeletal modeling and           #
# simulation. See http://opensim.stanford.edu and the NOTICE file         #
# for more information. OpenSim is developed at Stanford University       #
# and supported by the US National Institutes of Health (U54 GM072970,    #
# R24 HD065690) and by DARPA through the Warrior Web program.             #
#                                                                         #
# Copyright (c) 2005-2017 Stanford University and the Authors             #
# Author(s): Chris Dembia, Carmichael Ong                                 #
#                                                                         #
# Licensed under the Apache License, Version 2.0 (the "License");         #
# you may not use this file except in compliance with the License.        #
# You may obtain a copy of the License at                                 #
# http://www.apache.org/licenses/LICENSE-2.0.                             #
#                                                                         #
# Unless required by applicable law or agreed to in writing, software     #
# distributed under the License is distributed on an "AS IS" BASIS,       #
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or         #
# implied. See the License for the specific language governing            #
# permissions and limitations under the License.                          #
# ----------------------------------------------------------------------- #

# This script is intended to be used to combine entries for the 
# Tug Of War Competition.
# http://simtk-confluence.stanford.edu:8080/display/OpenSim/Designing+a+Muscle+for+a+Tug-of-War+Competition
# This script only works with Millard muscles.

import opensim as osim
import matplotlib.pyplot as plt
import re

import matplotlib
matplotlib.use('TkAgg')
import numpy as np

#from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
#from matplotlib.figure import Figure

#PYTHON GUI CODE SECTION
#--------------------

from tkinter import *

#ALL units are SI

newModel = osim.Model("C:/OpenSim 3.3/Models/Tug_of_War/Tug_of_War_Millard.osim")
newMuscleBase = newModel.getMuscles().get(0)
newMuscle = osim.Millard2012EquilibriumMuscle.safeDownCast(newMuscleBase)
initState = newModel.initSystem()

root = Tk()
root.title("paramter changer window")


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
global_figure_counter = 0

def alphaCom(num):
    newMuscle.set_pennation_angle_at_optimal(float(num))
    

def actCom(num):
    diff = deAct.get() - float(num)
    if diff < 0.03:
        setter = deAct.get() + (0.03 - diff)
        scaleDeAct.set(setter)
    if diff > 0.04:
        setter = deAct.get() - (diff - 0.04)
        scaleDeAct.set(setter)
    newMuscle.set_activation_time_constant(float(num))
    
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
    newModel.updConstraintSet()
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
	vector = osim.Vec3(0, 0.05, groundLoc)
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
    #print("tendon slack length: " + str(newMuscle.get_tendon_slack_length()))
    
#def file_len(fname):
#    count = newMuscle.set_tendon_slack_length(float(num))
#    blockLoc = newMuscle.getGeometryPath().getPathPointSet().get(1).getLocation().get(2)
#    # Left muscle
#    if blockLoc > 0:
#		groundLoc = blockLoc + fibLen + slackLen
#	# Right muscle
#    else:
#        groundLoc = blockLoc - fibLen - slackLen
#	vector = osim.Vec3(0, 0.05, groundLoc)
#	newMuscle.updGeometryPath().getPathPointSet().get(0).setLocation(initState, vector)
#    for line in fname:
#        count+=1
#    return count


Label(root, text = "Angle in degrees: ").pack()
scaleAlpha = Scale(root, variable = alpha, orient = HORIZONTAL, from_ = 0, to = 30, command = alphaCom)
scaleAlpha.pack()
scaleAlpha.set(0.)

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


labelTotLen = Label(root, text = "distance between origin and attachment: " + str(lMT))
labelTotLen.pack()

newFWDFile = open("_FWD_dynamics.xml",'w')

newFWDFile.write('<?xml version="1.0" encoding="UTF-8" ?>\n')
newFWDFile.write('<OpenSimDocument Version="30000">\n')
newFWDFile.write('\t\t<!--Name of the .osim file used to construct a model.-->\n')
newFWDFile.write('\t\t<model_file />\n')
newFWDFile.write('\t\t<!--Replace the model\'s force set with sets specified in <force_set_files>? If false, the force set is appended to.-->\n')
newFWDFile.write('\t\t<replace_force_set>false</replace_force_set>\n')
newFWDFile.write('\t\t<!--List of xml files used to construct an force set for the model.-->\n')
newFWDFile.write('\t\t<force_set_files />\n')
newFWDFile.write('\t\t<!--Directory used for writing results.-->\n')
newFWDFile.write('\t\t<results_directory>.</results_directory>\n')
newFWDFile.write('\t\t<!--Output precision.  It is 8 by default.-->\n')
newFWDFile.write('\t\t<output_precision>8</output_precision>\n')
newFWDFile.write('\t\t<!--Initial time for the simulation.-->\n')
newFWDFile.write('\t\t<initial_time>0</initial_time>\n')
newFWDFile.write('\t\t<!--Final time for the simulation.-->\n')
newFWDFile.write('\t\t<final_time>1</final_time>\n')
newFWDFile.write('\t\t<!--Flag indicating whether or not to compute equilibrium values for states other than the coordinates or speeds.  For example, equilibrium muscle fiber lengths or muscle forces.-->\n')
newFWDFile.write('\t\t<solve_for_equilibrium_for_auxiliary_states>true</solve_for_equilibrium_for_auxiliary_states>\n')
newFWDFile.write('\t\t<!--Maximum number of integrator steps.-->\n')
newFWDFile.write('\t\t<maximum_number_of_integrator_steps>20000</maximum_number_of_integrator_steps>\n')
newFWDFile.write('\t\t<!--Maximum integration step size.-->\n')
newFWDFile.write('\t\t<maximum_integrator_step_size>1</maximum_integrator_step_size>\n')
newFWDFile.write('\t\t<!--Minimum integration step size.-->\n')
newFWDFile.write('\t\t<minimum_integrator_step_size>1e-008</minimum_integrator_step_size>\n')
newFWDFile.write('\t\t<!--Integrator error tolerance. When the error is greater, the integrator step size is decreased.-->\n')
newFWDFile.write('\t\t<integrator_error_tolerance>1e-005</integrator_error_tolerance>\n')
newFWDFile.write('\t\t<!--Set of analyses to be run during the investigation.-->\n')
newFWDFile.write('\t\t<AnalysisSet name="Analyses">\n')
newFWDFile.write('\t\t\t<objects />\n')
newFWDFile.write('\t\t\t<groups />\n')
newFWDFile.write('\t\t</AnalysisSet>\n')
newFWDFile.write('\t\t<!--Controller objects in the model.-->\n')
newFWDFile.write('\t\t<ControllerSet name="Controllers">\n')
newFWDFile.write('\t\t\t<objects>\n')
newFWDFile.write('\t\t\t\t<ControlSetController>\n')
newFWDFile.write('\t\t\t\t\t<!--A Storage (.sto) or an XML control nodes file containing the controls for this controlSet.-->\n')
newFWDFile.write('\t\t\t\t</ControlSetController>\n')
newFWDFile.write('\t\t\t</objects>\n')
newFWDFile.write('\t\t\t<groups />\n')
newFWDFile.write('\t\t</ControllerSet>\n')
newFWDFile.write('\t\t<!--XML file (.xml) containing the forces applied to the model as ExternalLoads.-->\n')
newFWDFile.write('\t\t<external_loads_file />\n')
newFWDFile.write('\t\t<!--Storage file (.sto) containing the initial states for the forward simulation. This file often contains multiple rows of data, each row being a time-stamped array of states. The first column contains the time.  The rest of the columns contain the states in the order appropriate for the model. In a storage file, unlike a motion file (.mot), non-uniform time spacing is allowed.  If the user-specified initial time for a simulation does not correspond exactly to one of the time stamps in this file, inerpolation is NOT used because it is usually necessary to being a simulation from an exact set of states.  Instead, the closest earlier set of states is used. Having a states file that contains the entire trajectory of a simulations allows for corrective springs for perturbation analysis to be added.-->\n')
newFWDFile.write('\t\t<!--Flag (true or false) indicating whether or not the integrator should use a particular time stepping.  If true, the time stepping is extracted from the initial states file.  In this situation, therefore, the initial states file must contain all the time steps in a simulation and be written out to high precision (usually 20 decimal places).  Setting this flag to true can be useful when reproducing a previous forward simulation with as little drift as possible.  If this flag is false, the integrator is left to determine its own time stepping.-->\n')
newFWDFile.write('\t\t<use_specified_dt>false</use_specified_dt>\n')
newFWDFile.write('\t</ForwardTool>\n')
newFWDFile.write('</OpenSimDocument>\n')

newFWDFile.close()
    
def callBack(): 
    newModel.printToXML("C:/Users/mhmdk/Desktop/Co-op files/co-op semester 1/Python Code/newVersion.osim")
    newVersionModel = osim.Model("C:/Users/mhmdk/Desktop/Co-op files/co-op semester 1/Python Code/newVersion.osim")
    newVersionBase = newVersionModel.getMuscles().get(0)
    print(newVersionBase.getName())
    newVersionMuscle = osim.Millard2012EquilibriumMuscle.safeDownCast(newVersionBase)
    
    fibLen = float(lo.get())
    slackLen = float(ls.get())

    blockLoc = newVersionMuscle.getGeometryPath().getPathPointSet().get(1).getLocation().get(2)
    # Left muscle
    if blockLoc > 0:
        groundLoc = blockLoc + fibLen + slackLen
        print("fibLen: " + str(fibLen))
        print("slackLen: " + str(slackLen))
        print("blockLoc: " + str(blockLoc))
	# Right muscle
    else:
        groundLoc = blockLoc - fibLen - slackLen
	vector = osim.Vec3(0, 0.05, groundLoc)
    print(groundLoc)
    print("fibLen: " + str(fibLen))
    print("slackLen: " + str(slackLen))
    print("blockLoc: " + str(blockLoc))
    newVersionMuscle.updGeometryPath().getPathPointSet().get(0).setLocation(initState, vector)
    
    newVersionModel.initSystem()
    
    new_states_file = open("corrected_tug_of_war_states.sto", "w")

    new_states_file.write("Tug_of_War_Competition\n")
    new_states_file.write("nRows = 1\n")
    new_states_file.write("nColumns = 7\n")
    new_states_file.write("inDegree = no\n")
    new_states_file.write("endheader\n")
    new_states_file.write("time\tblock_tz\tblock_tz_u\tRightMuscle.activation\tRightMuscle.fiber_length\tLeftMuscle.activation\tLeftMuscle.fiber_length\n")
    new_states_file.write("0\t0\t0\t0.01\t" + str(newVersionMuscle.get_optimal_fiber_length()) +"\t0.01\t0.1")
    
    new_states_file.close()
    
    reporter = osim.ForceReporter(newVersionModel)
    newVersionModel.addAnalysis(reporter)
    fwd_tool = osim.ForwardTool()
    fwd_tool.setControlsFileName("C:\OpenSim 3.3\Models\Tug_of_War\Tug_of_War_Millard_controls_corrected.xml")
    fwd_tool.setStatesFileName("C:\Users\mhmdk\Desktop\Co-op files\co-op semester 1\Python Code\corrected_tug_of_war_states.sto")
    fwd_tool.setModel(newVersionModel)
    fwd_tool.run()
    updater()

button = Button(root, text = "GO!", command = callBack, font = "Verdana 9 bold", relief = RAISED, bg = "darkgreen", fg = "ivory")
button.pack()

timeLst = []
rightLst = []
leftLst = []
opener =  open("C:\Users\mhmdk\Desktop\Co-op files\co-op semester 1\Python Code\_ForceReporter_forces.sto")
count = 0
for line in opener:
    # line 15 is the line where the results start
    if count > 15: 
        temp = re.split("\t", line)
        rawTime = temp[0]
        rawRight = temp[1]
        rawLeft = temp[2]
        
        time = rawTime.strip("' ,\n\t")
        right = rawRight.strip("' ,\n\t")
        left = rawLeft.strip("' ,\n\t")
        
        timeLst.append(float(time))
        rightLst.append(float(right))
        leftLst.append(float(left))
        
    count+=1
opener.close()   

import matplotlib as mpl
import numpy as np
import sys
if sys.version_info[0] < 3:
    import Tkinter as tk
else:
    import tkinter as tk
import matplotlib.backends.tkagg as tkagg
from matplotlib.backends.backend_agg import FigureCanvasAgg


def draw_figure(canvas, figure, loc=(0, 0)):
    """ Draw a matplotlib figure onto a Tk canvas

    loc: location of top-left corner of figure on canvas in pixels.
    Inspired by matplotlib source: lib/matplotlib/backends/backend_tkagg.py
    """
    figure_canvas_agg = FigureCanvasAgg(figure)
    figure_canvas_agg.draw()
    figure_x, figure_y, figure_w, figure_h = figure.bbox.bounds
    figure_w, figure_h = int(figure_w), int(figure_h)
    photo = tk.PhotoImage(master=canvas, width=figure_w, height=figure_h)
    # Position: convert from top-left anchor to center anchor
    canvas.create_image(loc[0] + figure_w/2, loc[1] + figure_h/2, image=photo)

    # Unfortunately, there's no accessor for the pointer to the native renderer
    tkagg.blit(photo, figure_canvas_agg.get_renderer()._renderer, colormode=2)

    # Return a handle which contains a reference to the photo object
    # which must be kept live or else the picture disappears
    return photo

def updater():
    timeLst = []
    rightLst = []
    leftLst = []
    opener =  open("C:\Users\mhmdk\Desktop\Co-op files\co-op semester 1\Python Code\_ForceReporter_forces.sto")
    count = 0
    for line in opener:
        # line 15 is the line where the results start
        if count > 15: 
            temp = re.split("\t", line)
            rawTime = temp[0]
            rawRight = temp[1]
            rawLeft = temp[2]
            
            time = rawTime.strip("' ,\n\t")
            right = rawRight.strip("' ,\n\t")
            left = rawLeft.strip("' ,\n\t")
            
            timeLst.append(float(time))
            rightLst.append(float(right))
            leftLst.append(float(left))
            
        count += 1
    opener.close()
    
    # Create a canvas
    w, h = 500, 300
    window = Tk()
    window.title("plot")
    window.config(background='white')
    window.geometry("500x350")
    canvas = tk.Canvas(window, width=w, height=h, bg = "white")
    canvas.pack()
    
    # Generate some example data
    X = timeLst
    Y = rightLst
        
    # Create the figure we desire to add to an existing canvas
    fig = mpl.figure.Figure(figsize=(4, 3))
    ax = fig.add_axes([0, 0, 1, 1])
#    ax.set_axis_label("Force vs. time")
#    ax.set_axis_on()
    ax.set_xlabel("X axis")
    ax.set_ylabel("Y axis")
    ax.grid()
    ax.plot(X, Y, color = "red")
    
    # Keep this handle alive, or else figure will disappear
    fig_x, fig_y = 50, 50
    fig_photo = draw_figure(canvas, fig, loc=(fig_x, fig_y))
    fig_w, fig_h = fig_photo.width(), fig_photo.height()
    
    window.mainloop()
                
root.mainloop()