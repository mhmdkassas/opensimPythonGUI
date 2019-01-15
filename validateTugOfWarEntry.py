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

# This script is intended to be used to validate entries for the 
# Tug Of War Competition.
# http://simtk-confluence.stanford.edu:8080/display/OpenSim/Designing+a+Muscle+for+a+Tug-of-War+Competition
# This script only works with Millard muscles.

import math # for pi
import org.opensim.utils as utils

# LOAD MODEL, EXCITATIONS.
# ------------------------
# TODO maybe this should take a forward setup xml file?? Is that more difficult
# for students to use?

# TODO create dialog box (GUI) prompt to for a filename to a model, then load the model
# with loadModel().
model = getCurrentModel()
print(model.getName())

if not model:
    print("ERROR: Need to load a model first.")
else:
    print("Loaded model %s." % model.getName())

controlsFpath = utils.FileUtils.getInstance().browseForFilename(".xml",
    "Select the controls XML file.", 1)
controlSet = modeling.ControlSet(controlsFpath)
    
initState = model.initSystem()


# MUSCLE PROPERTIES
# -----------------
specificTension = 35 # N/cm^2
volumeMax = 100 # cm^3
powerMax = 175 # Watts
muscleTendonLengthMin = 0.15 # m
mucsleTendonLengthMax = 0.45 # m
tendonSlackLengthMin  = 0.10 # m
optimalFiberLengthMin = 0.05 # m
optimalFiberLengthMax = 0.20 # m
pennationAngleAtOptimalMin = 0
pennationAngleAtOptimalMax = 30 # degrees
maxContractionVelocityMin =  2 # optimal fiber lengths per second
maxContractionVelocityMax = 10 # optimal fiber lengths per second
excitationIntegralMax = 0.5
activationTimeConstantMin   = 0.010 # s
activationTimeConstantMax   = 0.020 # s
deactivationTimeConstantMin = 0.040 # s
deactivationTimeConstantMax = 0.060 # s
deactMinusActivTimeConstMin = 0.030 # s
deactMinusActivTimeConstMax = 0.040 # s

def check_inequality(min, val, max, name, units):
    if not (min <= val and val <= max):
        return "For %s, the inequailty %g %s <= %g %s <= %g %s does not hold.\n" % (
            name, min, units, val, units, max, units)
    else:
        return ""
        
def check(condition, errorMessage):
    if not condition:
        return "    " + errorMessage + "\n"
    else:
        return ""
        
def validateMuscle(baseMuscle):
    msg = '' # Store all error messages in this string.
    name = baseMuscle.getName()
    muscle = modeling.Millard2012EquilibriumMuscle.safeDownCast(baseMuscle)
    if not muscle:
        raise Exception("Muscle '%s' is not a Millard2012EquilibriumMuscle." % name)
    # print("Validating muscle '%s'." % name)
    area = muscle.get_max_isometric_force() / specificTension # in cm^2.
    
    # VOLUME
    # Fiber length should be expressed in centimeters.
    volume = area * (muscle.get_optimal_fiber_length() * 100)
    msg += check(volume <= volumeMax, "Volume of %g cm^3 exceeds maximum of %g cm^3." % (volume, volumeMax))
    
    # POWER
    power = (muscle.get_max_isometric_force() *
             muscle.get_max_contraction_velocity() *
             muscle.get_optimal_fiber_length()) # in Watts.
    msg += check(power <= powerMax, "Power of %g W exceeds maximum of %g W." % (power, powerMax))
    
    # TENDON SLACK LENGTH
    msg += check(muscle.get_tendon_slack_length() >= tendonSlackLengthMin, 
            "Tendon slack length of %g m is below minimum of %g m." % (
            muscle.get_tendon_slack_length(), tendonSlackLengthMin))
            
    # OPTIMAL FIBER LENGTH
    optFibLen = muscle.get_optimal_fiber_length()
    msg += check_inequality(optimalFiberLengthMin, optFibLen, optimalFiberLengthMax,
        "optimal fiber length", "m")    
        
    # PENNATION
    penn = muscle.get_pennation_angle_at_optimal() * 180.0 / math.pi # convert to degrees.
    msg += check_inequality(pennationAngleAtOptimalMin, penn, pennationAngleAtOptimalMax,
        "pennation angle at optimal fiber length", "deg")
        
    # MAX CONTRACTION VELOCITY
    maxVel = muscle.get_max_contraction_velocity()
    msg += check_inequality(maxContractionVelocityMin, maxVel, maxContractionVelocityMax,
        "max contraction velocity", "lM0/s")
  
    # ACTIVATION TIME CONSTANT
    activ = muscle.get_activation_time_constant()
    msg += check_inequality(activationTimeConstantMin, activ, activationTimeConstantMax,
        "activation time constant", "s")
    
    # DEACTIVATION TIME CONSTANT
    deact = muscle.get_deactivation_time_constant()
    msg += check_inequality(deactivationTimeConstantMin, deact, deactivationTimeConstantMax,
        "deactivation time constant", "s")
        
    # DEACTIVATION MINUS ACTIVATION TIME CONSTANT
    timeDiff = deact - activ
    msg += check_inequality(deactMinusActivTimeConstMin, timeDiff, deactMinusActivTimeConstMax,
        "[deact. - activ.] time constant", "s")
        
    # INITIAL MUSCLE-TENDON LENGTH
    initLength = muscle.getGeometryPath().getLength(initState)
    expectedInitLength = optFibLen + muscle.get_tendon_slack_length()
    msg += check(abs(initLength - expectedInitLength) < 1e-10,
        "Initial MTU length is %g m, but should be %g m + %g m = %g m." % (initLength,
        optFibLen, muscle.get_tendon_slack_length(), expectedInitLength))
    
    # TODO I don't know how crazy we want to get...should we check that there are only two path
    # points and that their x and y locations are correct?
    if msg != "":
        print("Muscle '%s' failed validation for the following reason(s):\n%sThe anti-doping agency has put you on their watch list.\n" % (name, msg))
    else:
        print("Muscle '%s' passed the validation!\n" % (name))

# Check each muscle in the model.
for imusc in range(model.getMuscles().getSize()):
    validateMuscle(model.getMuscles().get(imusc))


# EXCITATION SIGNAL
# -----------------

def integrate(x, y):
    assert(len(x) == len(y))
    integral = 0
    for i in range(len(x) - 1):
        width = x[i+1] - x[i]
        avgHeight = 0.5 * (y[i+1] + y[i])
        integral += width * avgHeight
    return integral

def clampExcitation(excitation):
    if excitation < 0:
        return 0
    elif excitation > 1:
        return 1
    else:
        return excitation

def validateControl(baseControl):
    msg = ''
    name = baseControl.getName()
    control = modeling.ControlLinear.safeDownCast(baseControl)
    if not control:
        raise Exception("Control %s must be a ControlLinear." % name)
    time = list()
    excitation = list()
    nodes = control.getControlValues()
    for inode in range(nodes.getSize()):
        time.append(nodes.get(inode).getTime())
        excitation.append(clampExcitation(nodes.get(inode).getValue()))
    # Check if time points are between 0 and 1
    msg += check((min(time) >= 0 and max(time) <= 1),
        "Excitation times for %s must be between 0 and 1" % (name))
    # Check if first and last time points are 0 and 1, respectively
    msg += check(abs(time[0]) < 1e-10, "First time point must be at 0s")
    msg += check(abs(time[-1] - 1.0) < 1e-10, "Last time point must be at 1s")
    # Check if time is sorted in increasing order
    msg += check(time == sorted(time),
        "Excitation times for %s must be in ascending order." % (name))
    # This is the exact integral, because we have a piecewise linear function.
    integral = integrate(time, excitation)
    msg += check(integral <= excitationIntegralMax, 
        "Excitation integral of %g exceeds maximum of %g." % (integral, excitationIntegralMax))
    
    if msg != "":
        print("Excitation '%s' failed validation:\n%s\n" % (name, msg))
    else:
        print("Excitation '%s' passed the validation!\n" % (name))
        
    
for ic in range(controlSet.getSize()):
    validateControl(controlSet.get(ic))
    








