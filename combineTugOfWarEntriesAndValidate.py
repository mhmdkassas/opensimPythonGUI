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

import math # for pi
import org.opensim.utils as utils
import string

# COMBINE OSIM MODELS
#--------------------
#
# Takes original millard model as the base and replaces properties in both
# models with properties taken from loaded user models.

baseModel = modeling.Model(getInstallDir()+"/Models/Tug_of_War/Tug_of_War_Millard.osim")
print(baseModel.getName())

if not baseModel:
    print("ERROR: Need to load base Tug_of_War_Millard.osim model first.")
else:
    print("Loaded model %s." % baseModel.getName())	

# Copy the base model to adjust
newModel = modeling.Model(baseModel)
newModel.initSystem();

# Load user right muscle
rightModelPath = utils.FileUtils.getInstance().browseForFilename(".osim",
    "Select the right user's Model file.", 1)
rightModel = modeling.Model(rightModelPath)
initState = rightModel.initSystem();

# Grab user for the right muscle (assumed naming convention is <user>_<file>.<type>)
rightModelFname = string.split( rightModelPath, '\\' )[-1];
rightModelUser = string.split( rightModelFname, '_' )[0];
print( rightModelUser )

# User should have muscle in the LEFT
userRightMuscleBase = rightModel.getMuscles().get('LeftMuscle')
userRightMuscle = modeling.Millard2012EquilibriumMuscle.safeDownCast(userRightMuscleBase)
print(userRightMuscle.getName())

# Get new right muscle to copy over parameters
rightMuscleBase = newModel.getMuscles().get(0)
rightMuscle = modeling.Millard2012EquilibriumMuscle.safeDownCast(rightMuscleBase)
print(rightMuscle.getName())

# Load user left muscle
leftModelPath = utils.FileUtils.getInstance().browseForFilename(".osim",
    "Select the left user's Model file.", 1)
leftModel = modeling.Model(leftModelPath)
leftModel.initSystem();

# Grab user for the left muscle (assumed naming convention is <user>_<file>.<type>)
leftModelFname = string.split( leftModelPath, '\\' )[-1];
leftModelUser = string.split( leftModelFname, '_' )[0];
print( leftModelUser )

# User should have muscle in the LEFT
userLeftMuscleBase = leftModel.getMuscles().get('LeftMuscle')
userLeftMuscle = modeling.Millard2012EquilibriumMuscle.safeDownCast(userLeftMuscleBase)
print(userLeftMuscle.getName())

# Get new left muscle to copy over parameters
leftMuscleBase = newModel.getMuscles().get(1)
leftMuscle = modeling.Millard2012EquilibriumMuscle.safeDownCast(leftMuscleBase)
print(leftMuscle.getName())

# Function for transfering properties...
def transferMuscleProperties( userMuscle, newMuscle, userName ):

	# Set new max isometric force
	newMuscle.set_max_isometric_force( userMuscle.get_max_isometric_force() )
	
	# Set new optimal fiber length
	newMuscle.set_optimal_fiber_length( userMuscle.get_optimal_fiber_length() )
	
	# Set new max contraction velocity
	newMuscle.set_max_contraction_velocity( userMuscle.get_max_contraction_velocity() )
	
	# Set new tendon slack length
	newMuscle.set_tendon_slack_length( userMuscle.get_tendon_slack_length() )
	
	# Set new pennation angle
	newMuscle.set_pennation_angle_at_optimal( userMuscle.get_pennation_angle_at_optimal() )
	
	# Set new activation time constant
	newMuscle.set_activation_time_constant( userMuscle.get_activation_time_constant() )
	
	# Set new deactivation time constant
	newMuscle.set_deactivation_time_constant( userMuscle.get_deactivation_time_constant() )
	
	# Set z-position of muscle on ground
	blockLoc = newMuscle.getGeometryPath().getPathPointSet().get(1).getLocation().get(2);
	
	# Left muscle
	if (blockLoc > 0):
		groundLoc = blockLoc + userMuscle.get_optimal_fiber_length() + userMuscle.get_tendon_slack_length()
	# Right muscle
	else:
		groundLoc = blockLoc - userMuscle.get_optimal_fiber_length() - userMuscle.get_tendon_slack_length()
	
	newMuscle.updGeometryPath().getPathPointSet().get(0).setLocation(initState, [0,0.05,groundLoc])
	
	# Name the muscle after the user
	newMuscle.setName(userName)

# Transfer properties!
transferMuscleProperties(userRightMuscle, rightMuscle, rightModelUser)
transferMuscleProperties(userLeftMuscle, leftMuscle, leftModelUser)

# Save model
modelName = string.join([rightModelUser, '_', leftModelUser, '_Tug_of_War_Millard.osim'], '');
newModel.print(modelName);

# COMBINE EXCITATIONS
#--------------------
#
# Takes excitation profiles from both users and combines them.
# Note, I gave up trying to do this in opensim and went with
# standard file io methods instead.

# Load right user control
rightControlFpath = utils.FileUtils.getInstance().browseForFilename(".xml",
    "Select the left user's controls XML file.", 1)

# Load left user control
leftControlFpath = utils.FileUtils.getInstance().browseForFilename(".xml",
    "Select the right user's controls XML file.", 1)

# Ensure same users
rightControlFname = string.split( rightControlFpath, '\\' )[-1];
rightControlUser = string.split( rightControlFname, '_' )[0];
assert( (rightControlUser in rightModelUser) and (rightModelUser in rightControlUser) )
print( rightControlUser )

leftControlFname = string.split( leftControlFpath, '\\' )[-1];
leftControlUser = string.split( leftControlFname, '_' )[0];
assert( (leftControlUser in leftModelUser) and (leftModelUser in leftControlUser) )
print( leftControlUser )

newControlFile = open(string.join([rightModelUser,'_',leftModelUser,'_Tug_of_War_Millard_control.xml'],''), 'w')

# Header
newControlFile.write('<?xml version="1.0" encoding="UTF-8" ?>\n')
newControlFile.write('<OpenSimDocument Version="30000">\n')
newControlFile.write('\t<ControlSet name="Control Set">\n')
newControlFile.write('\t\t<objects>\n')

# Read in muscle from each file
rightControlFile = open(rightControlFpath,'r')
ll = rightControlFile.readline()
while not ('LeftMuscle' in ll):
	ll = rightControlFile.readline()
	
newControlFile.write( string.join(['\t\t\t<ControlLinear name="',rightControlUser,'">\n'],'') )
ll = rightControlFile.readline()
while not ('</ControlLinear>' in ll):
	newControlFile.write( ll )
	ll = rightControlFile.readline()
newControlFile.write( ll )

leftControlFile = open(leftControlFpath,'r')
leftControlFile.readline()
while not ('LeftMuscle' in ll):
	ll = leftControlFile.readline()
	
newControlFile.write( string.join(['\t\t\t<ControlLinear name="',leftControlUser,'">\n'],'') )
ll = leftControlFile.readline()
while not ('</ControlLinear>' in ll):
	newControlFile.write( ll )
	ll = leftControlFile.readline()
newControlFile.write( ll )

newControlFile.write('\t\t</objects>\n')
newControlFile.write('\t\t<groups />\n')
newControlFile.write('\t</ControlSet>\n')
newControlFile.write('</OpenSimDocument>\n')

newControlFile.close()
leftControlFile.close()
rightControlFile.close()

# INITIAL STATES FILE
#--------------------
#
# Same thing, gave up trying to OpenSim it and just doing python file io

newStatesFile = open(string.join([rightModelUser,'_',leftModelUser,'_Tug_of_War_Millard_states.sto'],''),'w')

newStatesFile.write(string.join([rightModelUser,'_',leftModelUser,'_Tug_of_War_Millard_states\n'],''))
newStatesFile.write('version=1\n')
newStatesFile.write('nRows=1\n')
newStatesFile.write('nColumns=7\n')
newStatesFile.write('inDegrees=no\n')
newStatesFile.write('endheader\n')

# Header
newStatesFile.write('time\t')
newStatesFile.write('block_tz\t')
newStatesFile.write('block_tz_u\t')
newStatesFile.write(string.join([rightModelUser,'.activation\t'],''))
newStatesFile.write(string.join([rightModelUser,'.fiber_length\t'],''))
newStatesFile.write(string.join([leftModelUser,'.activation\t'],''))
newStatesFile.write(string.join([leftModelUser,'.fiber_length\n'],''))

# States
newStatesFile.write('0\t')
newStatesFile.write('0\t')
newStatesFile.write('0\t')
newStatesFile.write('0.01\t')
newStatesFile.write(string.join([str(rightMuscle.get_optimal_fiber_length()), '\t'],''))
newStatesFile.write('0.01\t')
newStatesFile.write(string.join([str(leftMuscle.get_optimal_fiber_length()), '\t'],''))

newStatesFile.close()

# FORWARD TOOL SETUP
#-------------------
#
# More io!

newFWDFile = open(string.join([rightModelUser,'_',leftModelUser,'_fwd_setup.xml'],''),'w')

newFWDFile.write('<?xml version="1.0" encoding="UTF-8" ?>\n')
newFWDFile.write('<OpenSimDocument Version="30000">\n')
newFWDFile.write(string.join(['\t<ForwardTool name="', rightModelUser, '_', leftModelUser, '">\n'],''))
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
newFWDFile.write(string.join(['\t\t\t\t\t<controls_file>',rightModelUser,'_',leftModelUser,'_Tug_of_War_Millard_control.xml</controls_file>\n'],''))
newFWDFile.write('\t\t\t\t</ControlSetController>\n')
newFWDFile.write('\t\t\t</objects>\n')
newFWDFile.write('\t\t\t<groups />\n')
newFWDFile.write('\t\t</ControllerSet>\n')
newFWDFile.write('\t\t<!--XML file (.xml) containing the forces applied to the model as ExternalLoads.-->\n')
newFWDFile.write('\t\t<external_loads_file />\n')
newFWDFile.write('\t\t<!--Storage file (.sto) containing the initial states for the forward simulation. This file often contains multiple rows of data, each row being a time-stamped array of states. The first column contains the time.  The rest of the columns contain the states in the order appropriate for the model. In a storage file, unlike a motion file (.mot), non-uniform time spacing is allowed.  If the user-specified initial time for a simulation does not correspond exactly to one of the time stamps in this file, inerpolation is NOT used because it is usually necessary to being a simulation from an exact set of states.  Instead, the closest earlier set of states is used. Having a states file that contains the entire trajectory of a simulations allows for corrective springs for perturbation analysis to be added.-->\n')
newFWDFile.write(string.join(['\t\t<states_file>',rightModelUser,'_',leftModelUser,'_Tug_of_War_Millard_states.sto</states_file>\n'],''))
newFWDFile.write('\t\t<!--Flag (true or false) indicating whether or not the integrator should use a particular time stepping.  If true, the time stepping is extracted from the initial states file.  In this situation, therefore, the initial states file must contain all the time steps in a simulation and be written out to high precision (usually 20 decimal places).  Setting this flag to true can be useful when reproducing a previous forward simulation with as little drift as possible.  If this flag is false, the integrator is left to determine its own time stepping.-->\n')
newFWDFile.write('\t\t<use_specified_dt>false</use_specified_dt>\n')
newFWDFile.write('\t</ForwardTool>\n')
newFWDFile.write('</OpenSimDocument>\n')

newFWDFile.close()

# RUN FORWARD SIMULATION
#-----------------------
#
# Back to OpenSim!
fwd_tool = modeling.ForwardTool(string.join([rightModelUser,'_',leftModelUser,'_fwd_setup.xml'],''))
fwd_tool.setModel(newModel)
fwd_tool.run()


# PLOT
# ----
loadModel(newModel)
states_fpath = string.join([rightModelUser, '_', leftModelUser, '_states.sto'], '')
loadMotion(states_fpath)


# CONTROLS
#plotterPanel = createPlotterPanel(rightModelUser + " vs " + leftModelUser + " controls")
# Load data from an external data source and plot
#src = addDataSource(plotterPanel, string.join([rightModelUser, '_', leftModelUser, '_controls.sto'], ''))
#crv2 = addCurve(plotterPanel, src, "time", rightModelUser)
#crv3 = addCurve(plotterPanel, src, "time", leftModelUser)
#crv2.setLegend(rightModelUser)
#crv3.setLegend(leftModelUser)


# POSITION
# Create a plotter panel and set the title
#plotterPanel = createPlotterPanel(rightModelUser + " (right, < 0) vs " + leftModelUser + " (left, > 0)")

# Load data from an external data source and plot
#src = addDataSource(plotterPanel, states_fpath)
#crv1 = addCurve(plotterPanel, src, "time", "block_tz")
#crv1.setLegend("block_tz")



# VALIDATE
# --------


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


controlSet = modeling.ControlSet(string.join([rightModelUser,'_',leftModelUser,'_Tug_of_War_Millard_control.xml'],''))
    
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
    


