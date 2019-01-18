# -*- coding: utf-8 -*-
"""
Created on Fri Jan 18 10:50:55 2019

@author: mhmdk
"""

import opensim as osim
from tkinter import*

newModel = osim.Model("C:/OpenSim 3.3/Models/Tug_of_War/Tug_of_War_Millard.osim")
newMuscleBase = newModel.getMuscles().get(0)
newMuscle = osim.Millard2012EquilibriumMuscle.safeDownCast(newMuscleBase)
initState = newModel.initSystem()

def main():
    root = Tk()
    
    Label(root, text = "Angle in degrees: ").pack()
    scaleAlpha = Scale(root, variable = alpha, orient = HORIZONTAL, from_ = 0, to = 30)
    scaleAlpha.pack()
    
    root.close()