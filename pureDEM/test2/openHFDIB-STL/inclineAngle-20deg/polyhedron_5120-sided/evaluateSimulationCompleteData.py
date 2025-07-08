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
    if not os.path.exists("postProc-results"):
        # if it doesn't exist, create it
        os.mkdir("postProc-results")
        
def exporting(path, header, time, data):
    """ template for exporting given data
    
        returns: data in text file
    """
    with open(path, "w") as file:
        file.write("\t".join(f"{h:<10}" for h in header) + "\n")

        columns = len(header) - 1

        for i in range(len(data)):
            row = [f"{float(time[i]):<10.6f}"]
            if columns == 1:
                row.append(f"{float(data[i]):<10.10f}")
            else:
                row.extend(f"{float(val):<10.10f}" for val in data[i][:columns])
            file.write("\t".join(row) + "\n")
        
class solidProperties:
    """ holds information about all bodies in the system
    
    """
    def __init__(self,bodyId,density,diameter,static = False):
        self.id = bodyId
        self.density = density
        self.diameter = diameter
        self.static = static
        self.particleParticleContactList = []
        self.particleWallContactList = []
        self.particlePositionList = []
        self.particleLinVelocityList = []
        self.particleAngularVelocityList = []
        self.inertia = []
        self.particleMassList = []
        self.particleMomentum = []
        self.particleAbsEnergy = []
        
        self.kinEnergy = []
        self.potEnergy = []
        self.mechEnergy = []
        self.angularEnergy = []
        self.disEnergy = []
        self.allContactEnergy = []
        
    def getRefenceStrings(self):
        return ["-- body %d ParticelMass  : "%(self.id),
                "-- body %d CoM                  : ("%(self.id),
                "-- body %d linear velocity      : ("%(self.id),
                "-- body %d angluar velocity     : "%(self.id)
                ]
    def getRefenceStringInertia(self):
        return "-- body %d moment of inerita magnitude : "%(self.id)
    
    def getActiveContactPair(self,secondBodyId,lastTime):
        for i in self.particleParticleContactList:
            if(i.tId == secondBodyId and i.isActive()):
                lastContactTime = i.getLatestContactTime()              
                if(lastContactTime != lastTime):
                    i.setActive(False)
                    return False
                else:
                    return i
        return False
    
    def getActiveContactWall(self,secondBodyId,lastTime):
        for i in self.particleWallContactList:
            if(i.wallId == secondBodyId and i.isActive()):
                lastContactTime = i.getLatestContactTime()
                if(lastContactTime != lastTime):
                    i.setActive(False)
                    return False
                else:
                    return i
        return False

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
            pe = self.particleMassList*np.dot(np.array([3.355217606024811, -9.21838460990976 ,0]),np.array([0.0,-0.26,0.0])-self.particlePositionList[i])
            self.potEnergy.append(pe)
    
    def calculateAllContactEnergy(self, TimeValues, solidsInDomain):
        fill_list = np.zeros(len(TimeValues[:-1]))
        for body in solidsInDomain:
            for j in body.particleParticleContactList:
                contact = 0
                if j.cId == self.id or j.tId == self.id:          
                    for i in range(len(TimeValues[:-1])):
                        if(j.contactTime[1] == TimeValues[i]):
                            j.ActiveReading = True 
                            j.contactTimeToPop = j.contactTime.copy()
                            j.contactTimeToPop.pop(0)
                            j.contactEnergyToPop = j.contactEnergy.copy()
                        if(j.ActiveReading):
                            if len(j.contactTimeToPop) > 0:
                                j.contactTimeToPop.pop(0)   
                                if solidsInDomain[j.cId].static == True or solidsInDomain[j.tId].static == True:
                                    contact = max(0,j.contactEnergyToPop[0])         
                                else:
                                    contact = max(0,j.contactEnergyToPop[0])/2   
                                j.contactEnergyToPop.pop(0)
                                contact += fill_list[i]
                                fill_list[i] = contact
        
        for j in self.particleWallContactList:
                contact = 0
                for i in range(len(TimeValues[:-1])):
                    if(j.contactTime[1] == TimeValues[i]):
                        j.ActiveReading = True
                        j.contactWallTimeToPop = j.contactTime.copy()
                        j.contactWallTimeToPop.pop(0)
                        j.contactWallEnergyToPop = j.contactWallEnergy.copy()
                        
                    if(j.ActiveReading):
                        if len(j.contactWallTimeToPop) > 0:
                            j.contactWallTimeToPop.pop(0)
                            contact = max(0,j.contactWallEnergyToPop[0])
                            j.contactWallEnergyToPop.pop(0)
                            contact += fill_list[i]
                            fill_list[i] = contact
                        else:
                            j.ActiveReading = False 
        self.allContactEnergy = fill_list

    def calculateMechanicEnergy(self,TimeValues, solidsInDomain):
        for i in range(len(self.kinEnergy)):
            e = self.kinEnergy[i] + self.potEnergy[i] + self.angularEnergy[i] + self.allContactEnergy[i]      
            self.mechEnergy.append(e)
    
    def calculateDisEnergy(self, TimeValues, solidsInDomain):
        fill_list = np.zeros(len(TimeValues[:-1]))
        for body in solidsInDomain:
            dis = 0  
            for j in body.particleParticleContactList:
                if j.cId == self.id or j.tId == self.id: 
                    for i in range(len(TimeValues[:-1])):
                        if(j.contactTime[1] == TimeValues[i]):
                            j.ActiveReading = True 
                            j.contactDisTimeToPop = j.contactTime.copy()
                            j.contactDisTimeToPop.pop(0)
                            j.contactDisEnergyToPop = j.contactDisEnergy.copy()
                        if(j.ActiveReading):
                            if len(j.contactDisTimeToPop) > 0:
                                j.contactDisTimeToPop.pop(0)
                                if solidsInDomain[j.cId].static == True or solidsInDomain[j.tId].static == True:
                                    dis = max(0,j.contactDisEnergyToPop[0])         
                                else:
                                    dis = max(0,j.contactDisEnergyToPop[0])/2
                                j.contactDisEnergyToPop.pop(0)
                                dis += fill_list[i]
                                if dis > 0:
                                    fill_list[i] = dis
                            else:
                                if dis > 0:
                                    fill_list[i] = dis

        for j in self.particleWallContactList:
            dis = 0
            for i in range(len(TimeValues[:-1])):
                if(j.contactTime[1] == TimeValues[i]):
                    j.ActiveReading = True 
                    j.contactWallTimeToPop = j.contactTime.copy()
                    j.contactWallTimeToPop.pop(0)
                    j.contactWallDisEnergyToPop = j.contactWallDisEnergy.copy()
                if(j.ActiveReading):
                    if len(j.contactWallTimeToPop) > 0:
                        j.contactWallTimeToPop.pop(0)
                        dis = max(0,j.contactWallDisEnergyToPop[0])         
                        j.contactWallDisEnergyToPop.pop(0)
                        dis += fill_list[i]
                        fill_list[i] = dis
                    else: 
                        fill_list[i] = dis
                    
        self.disEnergy = fill_list
            
    def calculateAbsEnergy(self):
        for i in range(len(self.mechEnergy)):
            absMech=self.mechEnergy[i]/self.mechEnergy[0]*100
            self.particleAbsEnergy.append(absMech)
            
    def Allexport(self, TimeValues, CFDtime):
        """ initializing data export

        """
        exporting("postProc-results/kin_energy_"+str(self.id), ["time", "kinetic_energy_body_"+ str(self.id)], TimeValues, self.kinEnergy)
        exporting("postProc-results/pot_energy_"+str(self.id), ["time", "potential_energy_body_"+ str(self.id)], TimeValues, self.potEnergy)
        exporting("postProc-results/mech_energy_"+str(self.id), ["time", "mechanical_energy_body_"+ str(self.id)], TimeValues, self.mechEnergy)
        exporting("postProc-results/dis_energy_"+str(self.id), ["time", "disipative_energy_body_"+ str(self.id)], TimeValues, self.disEnergy)
        exporting("postProc-results/abs_energy_"+str(self.id), ["time", "abs_energy_body_"+ str(self.id)], TimeValues, self.particleAbsEnergy)
        exporting("postProc-results/angular_energy_"+str(self.id), ["time", "angular_energy_body_"+ str(self.id)], TimeValues, self.angularEnergy)
        exporting("postProc-results/allContact_energy_"+str(self.id), ["time", "allContact_energy_body_"+ str(self.id)], TimeValues, self.allContactEnergy)
        exporting("postProc-results/momentum_"+str(self.id), ["time", "momentum_body_"+ str(self.id)], TimeValues, self.particleMomentum)
        exporting("postProc-results/position_body_"+str(self.id), ["time", "COM_x_body_"+ str(self.id), "COM_y_body_"+ str(self.id), "COM_z_body_"+ str(self.id)], TimeValues, self.particlePositionList)
        exporting("postProc-results/linVelocity_body_"+str(self.id), ["time", "linVel_x_body_"+ str(self.id), "linVel_y_body_"+ str(self.id), "linVel_z_body_"+ str(self.id)], TimeValues, self.particleLinVelocityList)
        exporting("postProc-results/angVelocity_body_"+str(self.id), ["time", "angVel_x_body_"+ str(self.id), "angVel_y_body_"+ str(self.id), "angVel_z_body_"+ str(self.id)], TimeValues, self.particleAngularVelocityList)
        exporting("postProc-results/inertia_body_"+str(self.id), ["time", "inertia_"+ str(self.id)], CFDtime, self.inertia)
    
class particleParticleContact:
    """ holds information about particle - particle contact of all bodies in the system
    
    """
    def __init__(self,contactPair):
        self.cId = contactPair[0] # body id of the fist body
        self.tId = contactPair[1] # body id of the second body
        self.contactTime = []
        self.contactTimeToPop = []
        self.contactForce = []
        self.contactNormal = []
        self.contactVolume = []
        self.contactArea = []
        self.delta = [] #overlap
        self.active = True
        self.diff_delta = []
        self.ActiveReading = False
        self.contactEnergyToPop = []
        self.contactEnergy = []
        
        self.contactDisTimeToPop = []
        self.contactDisEnergyToPop = []
        self.contactFNdclam = []
        self.contactDisEnergy = []

        self.fullContactDisEnergy = []

    def getLatestContactTime(self):
        return self.contactTime[-1]
    
    def isActive(self):
        return self.active
    
    def setActive(self,active):
        self.active = active
        
    def lastActive(self, TimeValues):
        if TimeValues[-2] == self.contactTime[-1]:
            self.active = True
        else: 
            self.active = False 

    def calculateContactEnergy(self):
        self.diff_delta = np.diff(self.delta,n=1)
        Ec = 0.0
        for i in range(len(self.diff_delta)):
            Ec +=(self.contactForce[i]+self.contactForce[i+1])*self.diff_delta[i]/2
            self.contactEnergy.append(Ec)  
            
    def calculateDisEnergy(self):
        self.diff_delta = np.diff(self.delta,n=1)
        Edis = 0.0
        for i in range(len(self.diff_delta)):
            Edis +=(self.contactFNdclam[i]+self.contactFNdclam[i+1])*abs(self.diff_delta[i])/2
            self.contactDisEnergy.append(Edis)  

    def adjust_dissipation_energy(self, TimeValues):
        dissipation_list = self.contactDisEnergy
        total_time_steps = len(TimeValues[:-1])
        for i in range(len(TimeValues[:-1])):
            if self.contactTime[1] == TimeValues[i]:
                initial_time_steps = i
        adjusted_list = [[0]] * initial_time_steps + [[value] for value in dissipation_list]
        if len(adjusted_list) < total_time_steps:
            last_value = dissipation_list[-1] if dissipation_list else 0
            adjusted_list += [[last_value]] * (total_time_steps - len(adjusted_list))
        else:
            adjusted_list = adjusted_list[:total_time_steps]
        self.fullContactDisEnergy = adjusted_list
             
    def Allexport(self):
        exporting("postProc-results/contact_dis_energy_"+str(self.cId)+"-"+str(self.tId), ["time", "contact_dis_energy_"+str(self.cId)+"-"+str(self.tId)], self.contactTime, self.contactDisEnergy)
        exporting("postProc-results/contact_energy_"+str(self.cId)+"-"+str(self.tId), ["time", "contact_energy"+str(self.cId)+"-"+str(self.tId)],  self.contactTime, self.contactEnergy)
        exporting("postProc-results/overlap_"+str(self.cId)+"-"+str(self.tId), ["time", "overlap_"+str(self.cId)+"-"+str(self.tId)],  self.contactTime, self.delta)
        exporting("postProc-results/contact_force_"+str(self.cId)+"-"+str(self.tId), ["time", "contact_force_"+str(self.cId)+"-"+str(self.tId)],  self.contactTime, self.delta)
        
class particleWallContact:
    """ holds information about wall - particle contact of all bodies in the system
    
    """
    def __init__(self,contactPair, point, normal):
        self.cId = contactPair[0]
        self.wallId = contactPair[1]
        self.contactTime = []
        self.contactWallTimeToPop = []
        self.contactForce = []
        self.contactNormal = []
        self.contactVolume = []
        self.contactArea = []
        self.delta = []
        self.active = True
        self.point = point
        self.normal = normal
        self.contactWallEnergy = []
        self.diff_delta=[]
        self.ActiveReading = False
        self.contactWallEnergyToPop = []

        self.contactWallFNdclam = []
        self.contactWallDisEnergy = []
        self.fullContactWallDisEnergy = []

    def getLatestContactTime(self):
        return self.contactTime[-1]
    
    def isActive(self):
        return self.active
    
    def setActive(self,active):
        self.active = active   

    def lastActive(self, TimeValues):
        if TimeValues[-2] == self.contactTime[-1]:
            self.active = True
        else: 
            self.active = False 

    def calculateWallEnergy(self):
        self.diff_delta = np.diff(self.delta,n=1)
        Ec = 0.0
        for i in range(len(self.diff_delta)):
            Ec +=(self.contactForce[i]+self.contactForce[i+1])*self.diff_delta[i]/2
            self.contactWallEnergy.append(Ec) 

    def calculateWallDisEnergy(self):
        self.diff_delta = np.diff(self.delta,n=1)
        Edis = 0.0
        for i in range(len(self.diff_delta)):
            Edis +=(self.contactWallFNdclam[i]+self.contactWallFNdclam[i+1])*abs(self.diff_delta[i])/2
            self.contactWallDisEnergy.append(Edis) 

    def adjust_dissipation_energy(self, TimeValues):
        dissipation_list = self.contactWallDisEnergy
        total_time_steps = len(TimeValues[:-1])
        for i in range(len(TimeValues[:-1])):
            if self.contactTime[1] == TimeValues[i]:
                initial_time_steps = i
        adjusted_list = [[0]] * initial_time_steps + [[value] for value in dissipation_list]
        if len(adjusted_list) < total_time_steps:
            last_value = dissipation_list[-1] if dissipation_list else 0
            adjusted_list += [[last_value]] * (total_time_steps - len(adjusted_list))
        else:
            adjusted_list = adjusted_list[:total_time_steps]
        self.fullContactWallDisEnergy = adjusted_list
            
    def Allexport(self):
        exporting("postProc-results/contact_dis_wallEnergy_"+str(self.cId), ["time", "contact_dis_wallEnergy_"+str(self.cId)], self.contactTime, self.contactWallDisEnergy)
        exporting("postProc-results/contact_wallEnergy_"+str(self.cId), ["time", "contact_wallEnergy_"+str(self.cId)], self.contactTime, self.contactWallEnergy)
        exporting("postProc-results/wall_overlap_"+str(self.cId), ["time", "wall_overlap_"+str(self.cId)],  self.contactTime, self.delta)
        exporting("postProc-results/wall_contactForce_"+str(self.cId), ["time", "wall_contactForce_"+str(self.cId)],  self.contactTime, self.contactForce)

def getWallDelta(solid):
    delta = 2*(solid.contactVolume[-1][0]/solid.contactArea[-1][0])
    return delta

def getParticleParticleDelta(solid):
    delta = 2*(solid.contactVolume[-1][0]/solid.contactArea[-1][0])
    return delta

def getPossibleContactPermutations(nSolids):
    possibleContactPermutations = []
    for i in range(nSolids):
        for j in range(i+1,nSolids):
            possibleContactPermutations.append([i,j])
    return possibleContactPermutations

def getPossibleWallContactPermutations(nSolids, nwall):
    possibleWallContactPermutations = []
    for j in range(1,nwall+1):
        for i in range(nSolids):
            possibleWallContactPermutations.append([i,j])
    return possibleWallContactPermutations

def countAmountContact(solidsInDomain):
    amount_PPContact = 0
    amount_PPContact_active = 0
    amount_WPContact = 0
    amount_WPContact_active = 0
    for body in solidsInDomain:
        amount_PPContact += len(body.particleParticleContactList)
        amount_WPContact += len(body.particleWallContactList)
        for contact in body.particleParticleContactList:
            if contact.active == True:
                amount_PPContact_active += 1
        for i in range(len(body.particleWallContactList)):
            if body.particleWallContactList[i].active == True:
                amount_WPContact_active += 1
            for j in range(i+1,len(body.particleWallContactList)):
                    contact_1 = body.particleWallContactList[i]
                    contact_2 = body.particleWallContactList[j]
                    if contact_1.contactTime[-1] == contact_2.contactTime[-1]:
                        amount_WPContact -= 1
                        if contact_1.active == True or contact_2.active == True:
                            amount_WPContact_active -= 1

    print("all partilce-particle contacts: %d, active: %d"%(amount_PPContact, amount_PPContact_active))
    print("all wall-particle contacts: %d, active: %d"%(amount_WPContact, amount_WPContact_active))
        
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
    """ extract body data from log.HFDIBDEMFoam for given reference string
    
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
                break

def getInertiaData(data,solid,initLine,endLine,time):  
    string = solid.getRefenceStringInertia()
    for i in range(initLine,endLine): 
        f_ind = data[i].find(string)
        if f_ind == 0:
            extracted_data = [float(strToVal) for strToVal in re.findall(r"-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *[-+]?\ *[0-9]+)?", data[i][f_ind+len(string)::])]
            solid.inertia.append(np.array(extracted_data)) 
            break            

def getParticleParticleContactData(data,contactPair,initLine,endLine,time,lastTime,solidsInDomain):
    """ extract particle-particle contact data from log.HFDIBDEMFoam for given reference string
    
        returns: data needed for evaluating particle-particle contact (delta, )
    
    """
    referenceString     = "-- Detected Particle-particle contact: -- body %d & -- body %d"%(contactPair[0],contactPair[1])
    stringList  = [
        "-- Particle-particle %d-%d contact FNe ("%(contactPair[0],contactPair[1]),
        "-- Particle-particle %d-%d contact normal ("%(contactPair[0],contactPair[1]),
        "-- Particle-particle %d-%d contact FNd clamped ("%(contactPair[0],contactPair[1]), 
        "-- Particle-particle %d-%d contact volume "%(contactPair[0],contactPair[1]),
        "-- Particle-particle %d-%d contact area "%(contactPair[0],contactPair[1]),
        ]
    for i in range(initLine,endLine):
        if referenceString in data[i] and data[i].split(referenceString)[1].strip() == "":
            contactPrt = solidsInDomain[contactPair[0]].getActiveContactPair(contactPair[1],lastTime)
            
            if(not contactPrt):
                solidsInDomain[contactPair[0]].particleParticleContactList.append(particleParticleContact(contactPair))
                contactPrt = solidsInDomain[contactPair[0]].particleParticleContactList[-1]

            contactPrt.contactTime.append(time)
            for j in range(i,endLine):
                for k in stringList:
                    f_ind = data[j].find(k)
                    if f_ind == 0:
                        extracted_data = [float(strToVal) for strToVal in re.findall(r"-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *[-+]?\ *[0-9]+)?", data[j][f_ind+len(k)::])]
                        if(k == stringList[0]):
                            contactPrt.contactForce.append(np.linalg.norm(np.array(extracted_data)))
                        elif(k == stringList[1]):
                            contactPrt.contactNormal = np.array(extracted_data)
                        elif (k == stringList[2]):
                            contactPrt.contactFNdclam.append(np.linalg.norm(np.array(extracted_data)))
                        elif (k == stringList[3]):
                            contactPrt.contactVolume.append(extracted_data)
                        elif (k == stringList[4]):
                            contactPrt.contactArea.append(extracted_data)
                        break
            contactPrt.delta.append(getParticleParticleDelta(contactPrt))
            break  

def getParticleWallContactData(data,contactPair,initLine,endLine,time,lastTime,solidsInDomain, wall):
    """ extract particle-wall contact data from log.HFDIBDEMFoam for given reference string
    
        returns: data needed for evaluating particle-wall contact (delta, contact force, contact normal)
    
    """
    referenceString = "-- Detected Particle-wall contact: -- body %d" % (contactPair[0]) 
    stringList = [
        "-- Particle-wall body %d contact FNe (" % (contactPair[0]),
        "-- Particle-wall body %d contact normal (" % (contactPair[0]),
        "-- Particle-wall body %d contact FNd clamped (" % (contactPair[0]),
        "-- Particle-wall body %d contact volume "%(contactPair[0]),
        "-- Particle-wall body %d contact area "%(contactPair[0]),
    ]
    
    for i in range(initLine, endLine):
        if referenceString in data[i] and data[i].split(referenceString)[1].strip() == "":
            contactPrt = solidsInDomain[contactPair[0]].getActiveContactWall(contactPair[1],lastTime)
            
            if(not contactPrt):
                solidsInDomain[contactPair[0]].particleWallContactList.append(particleWallContact(contactPair, wall.point, wall.normal))
                contactPrt = solidsInDomain[contactPair[0]].particleWallContactList[-1]
                
            contactPrt.contactTime.append(time)
            for j in range(i, endLine):
                for k in stringList:
                    f_ind = data[j].find(k)
                    if f_ind == 0:
                        extracted_data = [float(strToVal) for strToVal in re.findall(r"-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *[-+]?\ *[0-9]+)?", data[j][f_ind + len(k)::])]
                        if k == stringList[0]:
                            contactPrt.contactForce.append(np.linalg.norm(np.array(extracted_data)))
                        elif k == stringList[1]:
                            contactPrt.contactNormal = np.array(extracted_data)
                        elif (k == stringList[2]):
                            contactPrt.contactWallFNdclam.append(np.linalg.norm(np.array(extracted_data)))
                        elif (k == stringList[3]):
                            contactPrt.contactVolume.append(extracted_data)
                        elif (k == stringList[4]):
                            contactPrt.contactArea.append(extracted_data)
                        break

            if np.array_equal(contactPrt.normal, contactPrt.contactNormal):
                overlap = getWallDelta(contactPrt)
                contactPrt.delta.append(overlap)
            else:
                contactPrt.delta.append(0.0)

def main(filePath, maxNumberOfBodies):
    with open(filePath, 'r') as file:
        data = file.readlines()

    bodyNames, bodiesStatic = getBodyNames(data) #gets all bodies names and whether they are static or not
    nSolids, solidsInDomain = findBodiesInSimulation(data, maxNumberOfBodies, bodyNames, bodiesStatic)

    checkDir()
    contactPermutation = getPossibleContactPermutations(nSolids)
    nwall = 1
    contactWallPermutation = getPossibleWallContactPermutations(nSolids, nwall) #first number stands for body in contact with the wall, second number stands for wall id
    timeLines, TimeValues = getTimeStepLines(data)
    CFDStepLines, CFDtime = getCFDtime(data)
    for j in contactWallPermutation:
        wall = particleWallContact([j[0],1],np.array([0.0,0.0,0.0]),np.array([0.0,-1.0,0.0]))
 
        wallData=[wall]

    for i in range(len(timeLines[0:-1])): # to -1 because we exctract for future sequence (for integration)
        for j in solidsInDomain:
            getBodyData(data, j, timeLines[i], timeLines[i+1], TimeValues[i])
        for j in contactPermutation:
            getParticleParticleContactData(data, j, timeLines[i], timeLines[i+1], TimeValues[i], TimeValues[i-1], solidsInDomain)
        for j in contactWallPermutation:
            for wall in wallData:
                if j[1]==wall.wallId:
                    getParticleWallContactData(data, j, timeLines[i], timeLines[i+1], TimeValues[i], TimeValues[i-1], solidsInDomain, wall)
    
    #evaluate inertia data for each body
    for i in range(len(CFDStepLines[0:-1])):
        for j in solidsInDomain:
            getInertiaData(data, j, CFDStepLines[i], CFDStepLines[i+1], CFDtime[i])
    
    for body in solidsInDomain:
        #evaluate contact data for each body
        for contact in body.particleParticleContactList:
                contact.lastActive(TimeValues)
                contact.calculateContactEnergy()
                contact.calculateDisEnergy()
                contact.adjust_dissipation_energy(TimeValues)
                contact.Allexport()
        #evaluate wall contact data for each body
        for contact in body.particleWallContactList:
            if sum(contact.delta) != 0.0 and len(contact.contactForce) > 2:
                contact.lastActive(TimeValues)
                contact.calculateWallEnergy()
                contact.calculateWallDisEnergy()
                contact.adjust_dissipation_energy(TimeValues)
                contact.Allexport()
            else:
                contact.mark_for_deletion = True
        
        # remove contacts marked for deletion
        body.particleWallContactList = [contact for contact in body.particleWallContactList if not getattr(contact, 'mark_for_deletion', False)] 
        
    #count all contact in simulation
    # countAmountContact(solidsInDomain)

    #evaluate data for each body
    for body in solidsInDomain:
        if body.static == False:
            body.calculateKineticEnergy()
            body.calculatePotentialEnergy()
            body.calculateAllContactEnergy(TimeValues, solidsInDomain)
            body.calculateRotationEnergy(CFDtime, TimeValues)
            body.calculateMechanicEnergy(TimeValues, solidsInDomain)
            body.calculateMomentum()

            body.Allexport(TimeValues, CFDtime)
        else:
            body.calculateRotationEnergy(CFDtime, TimeValues)
            exporting("postProc-results/linVelocity_body_"+str(body.id), ["time", "linVel_x_body_"+ str(body.id), "linVel_y_body_"+ str(body.id), "linVel_z_body_"+ str(body.id)], TimeValues, body.particleLinVelocityList)
            exporting("postProc-results/position_body_"+str(body.id), ["time", "COM_x_body_"+ str(body.id), "COM_y_body_"+ str(body.id), "COM_z_body_"+ str(body.id)], TimeValues, body.particlePositionList)
            exporting("postProc-results/angular_energy_"+str(body.id), ["time", "angular_energy_body_"+ str(body.id)], TimeValues, body.angularEnergy)
            exporting("postProc-results/inertia_body_"+str(body.id), ["time", "inertia_"+ str(body.id)], CFDtime, body.inertia)

    systemEnergy = 0
    mechEnergy = 0
    AllDisEnergy = 0
    AllKinEnergy = 0
    AllPotEnergy = 0
    AllContactEnergy = 0
    AllAngularEnergy = 0
    for body in solidsInDomain:
        if len(body.mechEnergy) != 0:
            mechEnergy += np.array(body.mechEnergy)   
        if len(body.kinEnergy) != 0:
            AllKinEnergy += np.array(body.kinEnergy) 
        if len(body.potEnergy) != 0:
            AllPotEnergy += np.array(body.potEnergy)
        if len(body.allContactEnergy) != 0:
            AllContactEnergy += np.array(body.allContactEnergy)
        if len(body.angularEnergy) != 0:
            AllAngularEnergy += np.array(body.angularEnergy)
        for i in body.particleParticleContactList:
            if body.id == i.cId or body.id == i.tId:
                AllDisEnergy += np.array(i.fullContactDisEnergy)
        for i in body.particleWallContactList:
            AllDisEnergy += np.array(i.fullContactWallDisEnergy)

    systemEnergy = np.array(mechEnergy) + np.array(AllDisEnergy)
    energy = np.array(AllContactEnergy)+np.array(AllKinEnergy)

    exporting("postProc-results/all_mechanical_energy", ["time", "mechanical_energy"], TimeValues, systemEnergy)
    exporting("postProc-results/all_kinetic_energy", ["time", "kinetic_energy"], TimeValues, AllKinEnergy)
    exporting("postProc-results/all_potential_energy", ["time", "potential_energy"], TimeValues, AllPotEnergy)
    exporting("postProc-results/all_contact_energy", ["time", "contact_energy"], TimeValues, AllContactEnergy)
    exporting("postProc-results/all_angular_energy", ["time", "angular_energy"], TimeValues, AllAngularEnergy)
    exporting("postProc-results/all_kinContact_energy", ["time", "kinContact_energy"], TimeValues, energy)

data = main("log.HFDIBDEMFoam",35)