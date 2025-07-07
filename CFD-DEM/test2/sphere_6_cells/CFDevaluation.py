#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#=FILE DESCRIPTION======================================================
#
#Python script to extract data from HFDIBDEMFoam for each body in given time and calculates energy of the system and energy in contact.
#
#=======================================================================

#import libraries
import matplotlib.pyplot as plt
import numpy as np
import re
import os
import shutil

if os.path.exists("CFDpostProc-results"):
    shutil.rmtree("CFDpostProc-results")

def getBodyNames(data):
    """ finds all bodies in the system and assigns properties (static/not static)
        
        returns: list of body names and their property in the system

    """
    bodyNames = []
    referenceString = "Time = "
    soughtLine = 0
    for i in range(len(data)):
        f_ind = data[i].find(referenceString)
        if f_ind == 0:
            soughtLine = i
            break
    # -- find the number of bodies
    InitString = "Found bodyGeom for "
    endString = ", the body is: "

    for i in range(0,soughtLine):
        f_ind = data[i].find(InitString)
        if(f_ind == 0):          
            bodyNames.append(data[i][f_ind+len(InitString):data[i].find(endString)])
        if(data[i].find("Creating immersed body based on:") != -1):
            break

    bodiesStatic = bodyNames.copy()
    # print("bodyNames: ",bodyNames)
    for i in range(0,soughtLine):
        f_ind = data[i].find(" is a static body")
        if(f_ind != -1):
            name = data[i][:f_ind]
            bodiesStatic[bodiesStatic.index(name)] = True
            

    for i in bodiesStatic:
        if(type(i) == str):
            bodiesStatic[bodiesStatic.index(i)] = False
    # print("bodiesStatic: ",bodiesStatic)
    return bodyNames,bodiesStatic

def findBodiesInSimulation(data,maxLimit,bodyNames,bodiesStatic):
    """ calculate the number of bodies in the simulation
    
        returns: n of solids in simulation and list of solid properties (id, density, diameter)
    """
    solidsInDomain = []
    # -- estimate size of the simulation initialisation
    referenceString = "Time = "
    soughtLine = 0
    soughtLine2 = 0
    for i in range(len(data)):
        f_ind = data[i].find(referenceString)
        if f_ind == 0:
            soughtLine = i
            break

    for k in range(soughtLine,len(data)):
        f_ind = data[k].find(" DEM - CFD Time: ")
        if f_ind == 0:
            soughtLine2 = k
            # print(k)
            break
    # -- find the number of bodies    
    nFound = 0  
#    print ("soughtLine: ",soughtLine)
    for i in range(0,soughtLine):
        f_ind = data[i].find("New bodyID: ")
#        print(" i: ",i," f_ind: ",f_ind)
        if f_ind == 0:
            for j in range(nFound,maxLimit):
                if(data[i].find(getNewBodyString(j)) == 0):
#                    print(getNewBodyString(j))
                    nFound += 1
                    isStatic = False
                    for pos in range(len(bodyNames)):
                        if(data[i].find(bodyNames[pos]) != -1):
                            isStatic = bodiesStatic[pos]
                            break
                    extracted_data = [float(strToVal) for strToVal in re.findall(r"-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *[-+]?\ *[0-9]+)?", data[i][data[i].find("rhoS: ")+len("rhoS:  rho [1 -3 0 0 0 0 0]")::])]
                    
                    solidsInDomain.append(solidProperties(j,extracted_data[0],extracted_data[1],isStatic))
                    break
    #print("soughtLine: ",soughtLine," soughtLine2: ",soughtLine2," nFound: ",nFound)
    for i in range(soughtLine,soughtLine2):
        #print("i: ",i," data[i]: ",data[i])
        if(data[i].find(" current center of mass position: (") != -1):
            #print(data[i])
            for j in solidsInDomain: #j represents the body
                if(j.static and data[i].find(getStaticBodyString(j.id)) == 0):
                    
                    extracted_data = [float(strToVal) for strToVal in re.findall(r"-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *[-+]?\ *[0-9]+)?", data[i][data[i].find(getStaticBodyString(j.id))+len(getStaticBodyString(j.id))::])]
                    # appends the correcrt data to each properties in solidProperties
                    j.particlePositionList.append(np.array(extracted_data))
                    j.particleLinVelocityList.append(np.array([0.0,0.0,0.0]))
                    j.particleAngularVelocityList.append(np.array([0.0]))
                    break
                if (j.static and data[i].find(getStaticInertia(j.id)) == 0):
                    extracted_data1 = [float(strToVal) for strToVal in re.findall(r"-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *[-+]?\ *[0-9]+)?", data[i][data[i].find(getStaticInertia(j.id))+len(getStaticInertia(j.id))::])]
                    j.inertia.append(np.array(extracted_data1))
                    break

    return nFound,solidsInDomain
    # -- find the number of bodies

def getNewBodyString(bodyId):
    return "New bodyID: %d name: "%(bodyId)

def getStaticBodyString(bodyId):
    return "-- body: %d current center of mass position: ("%(bodyId)

def getStaticInertia(bodyId):
    return "-- body %d moment of inerita magnitude : "%(bodyId)

def checkDir():
    """ checks if postProc-result path exist, if not, creates one

    """
    if not os.path.exists("CFDpostProc-results"):
        # if it doesn't exist, create it
        os.mkdir("CFDpostProc-results")

def exporting(path, header, time, data):
    """ template for exporting given data
    
        returns: data in a text file
    """
    file_exists = os.path.isfile(path)  

    with open(path, "a") as file: 
        if not file_exists or os.stat(path).st_size == 0: 
            file.write("\t".join(f"{h:<10}" for h in header) + "\n")

        columns = len(header) - 1

        for i in range(len(data)):
            row = [f"{float(time[i]):<10.6f}"]

            if isinstance(data[i], (list, np.ndarray)):
                if len(data[i]) == 1:
                    row.append(f"{float(data[i][0]):<10.10f}") 
                else:
                    row.extend(f"{float(val):<10.10f}" for val in data[i][:columns])  
            else:
                row.append(f"{float(data[i]):<10.10f}")  
            file.write("\t".join(row) + "\n")

class solidProperties:
    """ holds information about all bodies in the system
    
    """
    def __init__(self,bodyId,density,diameter,static = False):
        self.id = bodyId
        self.density = density
        self.diameter = diameter
        self.static = static
        self.particlePositionList = []
        self.particleLinVelocityList = []
        self.particleAngularVelocityList = []
        self.particleFCouplingList = []
        self.particleFfluidList = []
        self.particleTransIncrementList = []
        self.inertia = []
        self.particleMassList = []
        self.particleMomentum = []
        self.diff_transIncrement = []
        self.Reynolds = []
        
        self.kinEnergy = []
        self.potEnergy = []
        self.mechEnergy = []
        self.angularEnergy = []
        self.dragEnergy = []
        self.fluidEnergy = []
        
    def getRefenceStrings(self):
        return ["-- body %d ParticelMass  : "%(self.id),
                "-- body %d CoM                  : ("%(self.id),
                "-- body %d linear velocity      : ("%(self.id),
                "-- body %d angluar velocity     : "%(self.id),
                "-- body %d translation increment: ("%(self.id),
                "-- body %d FCoupling: ("%(self.id),
                "-- body %d Ffluid: ("%(self.id),

                ]
    def getRefenceStringInertia(self):
        return "-- body %d moment of inerita magnitude : "%(self.id)

    def calculateRotationEnergy(self, CFDtime, TimeValues):
        time_intervals = int(len(TimeValues)/len(CFDtime))  
        for i in range(len(self.inertia)):
            inertia_value = self.inertia[i][0]
            for j in range(i*time_intervals, min(i*time_intervals + time_intervals, len(self.particleAngularVelocityList))):
                energy = 0.5*inertia_value*self.particleAngularVelocityList[j]**2
                self.angularEnergy.append(energy)

    def calculateMomentum(self):
        for i in range(len(self.particleLinVelocityList)):
            moment = self.particleMassList*np.linalg.norm(self.particleLinVelocityList[i])
            self.particleMomentum.append(moment)
    
    def calculateKineticEnergy(self):
        for i in range(len(self.particleLinVelocityList)):
            ke = 0.5 * self.particleMassList * np.linalg.norm(self.particleLinVelocityList[i]) ** 2
            self.kinEnergy.append(ke)
    
    def calculatePotentialEnergy(self):
        for i in range(len(self.particleLinVelocityList)):
            pe = self.particleMassList*np.dot(np.array([0,-9.81,0]),np.array([0.0,0.0,0.0])-self.particlePositionList[i])
            # pe = 0
            self.potEnergy.append(pe)

    def calculateMechanicEnergy(self):
        for i in range(len(self.dragEnergy)):
            e = self.angularEnergy[i] + self.dragEnergy[i] + self.fluidEnergy[i] + self.kinEnergy[i] + self.potEnergy[i]
            self.mechEnergy.append(e)

    def calculateDragEnergy(self):
        self.diff_transIncrement = np.diff(self.particleTransIncrementList,n=1)
        drag = 0.0
        for i in range(len(self.particleTransIncrementList)-1):
            drag +=(self.particleFCouplingList[i]+self.particleFCouplingList[i+1])*abs(self.particleTransIncrementList[i])/2
            self.dragEnergy.append(drag)

    def calculateFluidEnergy(self):
        self.diff_transIncrement = np.diff(self.particleTransIncrementList,n=1)
        fluid = 0.0
        for i in range(len(self.particleTransIncrementList)-1):
            fluid +=(self.particleFfluidList[i]+self.particleFfluidList[i+1])*abs(self.particleTransIncrementList[i])/2
            self.fluidEnergy.append(fluid)

    def calculateReynolds(self):
        nu = 1e-6
        for i in range(len(self.particleLinVelocityList)):
            re = self.diameter*np.linalg.norm(self.particleLinVelocityList[i])/nu
            self.Reynolds.append(re)
            
    def Allexport(self, TimeValues, CFDtime):
        """ initializing data export

        """
        print(self.particleMassList)
        exporting("CFDpostProc-results/kin_energy_"+str(self.id), ["time", "kinetic_energy_body_"+ str(self.id)], TimeValues, self.kinEnergy)
        exporting("CFDpostProc-results/pot_energy_"+str(self.id), ["time", "potential_energy_body_"+ str(self.id)], TimeValues, self.potEnergy)
        exporting("CFDpostProc-results/mech_energy_"+str(self.id), ["time", "mechanical_energy_body_"+ str(self.id)], TimeValues, self.mechEnergy)
        exporting("CFDpostProc-results/drag_energy_"+str(self.id), ["time", "drag_energy_body_"+ str(self.id)], TimeValues, self.dragEnergy)
        exporting("CFDpostProc-results/angular_energy_"+str(self.id), ["time", "angular_energy_body_"+ str(self.id)], TimeValues, self.angularEnergy)
        exporting("CFDpostProc-results/momentum_"+str(self.id), ["time", "momentum_body_"+ str(self.id)], TimeValues, self.particleMomentum)
        exporting("CFDpostProc-results/position_body_"+str(self.id), ["time", "COM_x_body_"+ str(self.id), "COM_y_body_"+ str(self.id), "COM_z_body_"+ str(self.id)], TimeValues, self.particlePositionList)
        exporting("CFDpostProc-results/linVelocity_body_"+str(self.id), ["time", "linVel_x_body_"+ str(self.id), "linVel_y_body_"+ str(self.id), "linVel_z_body_"+ str(self.id)], TimeValues, self.particleLinVelocityList)
        exporting("CFDpostProc-results/angVelocity_body_"+str(self.id), ["time", "angVel_x_body_"+ str(self.id), "angVel_y_body_"+ str(self.id), "angVel_z_body_"+ str(self.id)], TimeValues, self.particleAngularVelocityList)
        exporting("CFDpostProc-results/inertia_body_"+str(self.id), ["time", "inertia_"+ str(self.id)], CFDtime, self.inertia)
        exporting("CFDpostProc-results/FCoupling_body_"+str(self.id), ["time", "FCoupling_body_"+ str(self.id)], TimeValues, self.particleFCouplingList)
        exporting("CFDpostProc-results/Ffluid_body_"+str(self.id), ["time", "Ffluid_body_"+ str(self.id)], TimeValues, self.particleFfluidList)
        exporting("CFDpostProc-results/transIncrement_body_"+str(self.id), ["time", "translation_increment_body_"+ str(self.id)], TimeValues, self.particleTransIncrementList)
        exporting("CFDpostProc-results/reynolds_body_"+str(self.id), ["time", "reynolds_body_"+ str(self.id)], TimeValues, self.Reynolds)
        exporting("CFDpostProc-results/fluid_energy_"+str(self.id), ["time", "fluid_energy_body_"+ str(self.id)], TimeValues, self.fluidEnergy)
        
def getCFDtime(data):
    CFDStepLines = []
    CFDtime = []
    referenceString = "Time = "
    for i in range(len(data)):
        f_ind = data[i].find(referenceString)
        if f_ind == 0:
            CFDtime.append(float(data[i][f_ind+len(referenceString)::]))
            CFDStepLines.append(i)
    lenght = len(data)
    CFDStepLines.append(lenght)
    return CFDStepLines,CFDtime

def getTimeStepLines(data):
    """ finds where and what value has reference string

        returns: lines where reference string was found and all time values 
    """
    timeStepLines = []
    timeValues = []
    # -- estimate size of the simulation initialisation
    referenceString = " DEM - CFD Time: "
    for i in range(len(data)):
        f_ind = data[i].find(referenceString)
        if f_ind == 0:
            timeValues.append(float(data[i][f_ind+len(referenceString)::]))
            timeStepLines.append(i)

    return timeStepLines,timeValues

def getBodyData(data,solid,initLine,endLine,time):   #solid in solidsInDomain
    """ extract body data from log.pimpleHFDIBFoam for given reference string
    
        returns: data for each particle (velocity, position, mass)
    
    """
    stringList = solid.getRefenceStrings()
    for i in range(initLine,endLine):
        for j in stringList: 
            f_ind = data[i].find(j)
            if f_ind == 0:
                extracted_data = [float(strToVal) for strToVal in re.findall(r"-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *[-+]?\ *[0-9]+)?", data[i][f_ind+len(j)::])]
                if(j == stringList[0]):
                    solid.particleMassList = extracted_data[0]    
                elif(j == stringList[1]):
                    solid.particlePositionList.append(np.array(extracted_data))
                elif(j == stringList[2]):
                    solid.particleLinVelocityList.append(np.array(extracted_data))
                elif(j == stringList[3]):
                    solid.particleAngularVelocityList.append(np.array(extracted_data))
                elif(j == stringList[5]):
                    solid.particleFCouplingList.append(np.linalg.norm(np.array(extracted_data)))
                elif(j == stringList[4]):
                    solid.particleTransIncrementList.append(np.linalg.norm(np.array(extracted_data)))
                elif(j == stringList[6]):
                    solid.particleFfluidList.append(np.linalg.norm(np.array(extracted_data)))
                break

def getInertiaData(data,solid,initLine,endLine,time):  
    string = solid.getRefenceStringInertia()
    for i in range(initLine,endLine): 
        f_ind = data[i].find(string)
        if f_ind == 0:
            extracted_data = [float(strToVal) for strToVal in re.findall(r"-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *[-+]?\ *[0-9]+)?", data[i][f_ind+len(string)::])]
            solid.inertia.append(np.array(extracted_data)) 
            break            

def main(filePath, maxNumberOfBodies):
    with open(filePath, 'r') as file:
        data = file.readlines()

    bodyNames, bodiesStatic = getBodyNames(data) #gets all bodies names and whether they are static or not
    nSolids, solidsInDomain = findBodiesInSimulation(data, maxNumberOfBodies, bodyNames, bodiesStatic)

    checkDir()
    timeLines, TimeValues = getTimeStepLines(data)
    CFDStepLines, CFDtime = getCFDtime(data)

    for i in range(len(timeLines[0:-1])): # to -1 because we exctract for future sequence (for integration)
        for j in solidsInDomain:
            getBodyData(data, j, timeLines[i], timeLines[i+1], TimeValues[i])

    #evaluate inertia data for each body
    for i in range(len(CFDStepLines[0:-1])):
        for j in solidsInDomain:
            getInertiaData(data, j, CFDStepLines[i], CFDStepLines[i+1], CFDtime[i])

    #evaluate data for each body
    for body in solidsInDomain:
        if body.static == False:
            body.calculateKineticEnergy()
            body.calculatePotentialEnergy()
            body.calculateDragEnergy()
            body.calculateFluidEnergy()
            body.calculateRotationEnergy(CFDtime, TimeValues)
            body.calculateMechanicEnergy()
            body.calculateMomentum()
            body.calculateReynolds()

            body.Allexport(TimeValues, CFDtime)
        else:
            body.calculateRotationEnergy(CFDtime, TimeValues)
            exporting("CFDpostProc-results/linVelocity_body_"+str(body.id), ["time", "linVel_x_body_"+ str(body.id), "linVel_y_body_"+ str(body.id), "linVel_z_body_"+ str(body.id)], TimeValues, body.particleLinVelocityList)
            exporting("CFDpostProc-results/position_body_"+str(body.id), ["time", "COM_x_body_"+ str(body.id), "COM_y_body_"+ str(body.id), "COM_z_body_"+ str(body.id)], TimeValues, body.particlePositionList)
            exporting("CFDpostProc-results/angular_energy_"+str(body.id), ["time", "angular_energy_body_"+ str(body.id)], TimeValues, body.angularEnergy)
            exporting("CFDpostProc-results/inertia_body_"+str(body.id), ["time", "inertia_"+ str(body.id)], CFDtime, body.inertia)      

main("log.pimpleHFDIBFoam",5)