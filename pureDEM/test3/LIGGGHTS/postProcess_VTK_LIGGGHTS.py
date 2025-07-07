#!/usr/bin/python3

import re
import numpy as np
#import matplotlib.pyplot as plt
import os

def canBeConvertedToFloat(input):
    try:
        float(input)
        return True
    except ValueError:
        return False

def getTimeLevelList(caseDir,fileName):
    "Returns a list with the particle numbers in the case directory."
    Directory = caseDir       
    Strings   = list(set([numStr for numStr in os.listdir(Directory)]))

    Full_List = []
    l = len(fileName)
    for i in range(len(Strings)):
        if(Strings[i][-4:] == '.vtk' and Strings[i][:l]==fileName and canBeConvertedToFloat(Strings[i][l:-4])):
            h = Strings[i][:-4]
            Full_List.append(int(h[l:]))
    if(len(Full_List) > 0):            
        Full_List.sort(key = float)
    return Full_List 

class simulationResults:
    def __init__(self):
        "Class to store the data of a sphere body."
        self.timeLevels = {}
        self.globalIDList = {}

    def addTimeLevel(self, timeLevel, time):
        self.timeLevels[time] = timeLevel

    def addBodiesToGlobalList(self):
        timeList = self.getTimes()
        timeList.sort(key = float)
        latestTime = timeList[-1]
        for id in self.timeLevels[latestTime].getParticleList():
            if(id not in self.globalIDList):
                self.globalIDList[id] = id

    def getTimeLevel(self, time):
        try:
            return self.timeLevels[time]
        except:
            print("Time level not found")
    
    def getTimes(self):
        return list(self.timeLevels.keys())
    
    def getFirstTimeLevel(self):
        timeList = self.getTimes()
        timeList.sort(key = float)
        return self.timeLevels[timeList[0]]

    def getLatestTimeLevel(self):
        timeList = self.getTimes()
        timeList.sort(key = float)
        return self.timeLevels[timeList[-1]]

    def getBodyResults(self,id):
        timeList = self.getTimes()
        timeList.sort(key = float)
        
        if(not os.path.exists("processedResults")):
            os.system("mkdir processedResults")
        
        with open("processedResults/particle_"+str(id)+".dat","w") as file:
            headers = ["Time", "x","y","z", "vx", "vy", "vz", "fx" ,"fy", "fz", "omega", "radius", "mass" ]
            header_line = "{:<10}".format(headers[0]) + "\t" + "{:<20}".format(headers[1]) + "\t"
            header_line += "\t".join("{:<20}".format(header) for header in headers[2:]) + "\n"
            file.write(header_line)
            for key in timeList:
                if(not self.timeLevels[key].isParticlePresent(id)):
                    print("-- warning : Particle ", str(id),"was not found at time: ", key)
                    continue
                body = self.timeLevels[key].getParticle(id)
                body.calculateKineticEnergy()
                row_values =  ["{:<10}".format(str(key))] 
                row_values += ["{:<20.6f}".format(float(body.pos[i])) for i in range(3)]
                row_values += ["{:<20.6f}".format(float(body.linVel[i])) for i in range(3)]
                row_values += ["{:<20.6f}".format(float(body.force[i])) for i in range(3)]
                row_values += ["{:<20.6f}".format(float(np.linalg.norm(body.omega)))]
                row_values += ["{:<20.6f}".format(float(body.r))]
                row_values += ["{:<20.6f}".format(float(body.mass))]
                file.write("\t".join(row_values) + "\n")

        print("Results for particle ", str(id), "were written to processedResults/particle_"+str(id)+".dat")

    def getBodyNonZeroForces(self,id):
        timeList = self.getTimes()
        timeList.sort(key = float)
        
        if(not os.path.exists("processedResults")):
            os.system("mkdir processedResults")
        
        with open("processedResults/particleForces_"+str(id)+"_.dat","w") as file:
            headers = ["Time", "x","y","z", "vx", "vy", "vz", "fx" ,"fy", "fz", "omega", "radius", "mass" ,"kinetic"]
            header_line = "{:<10}".format(headers[0]) + "\t" + "{:<20}".format(headers[1]) + "\t"
            header_line += "\t".join("{:<20}".format(header) for header in headers[2:]) + "\n"
            
            file.write(header_line)
            for key in timeList:
                if(not self.timeLevels[key].isParticlePresent(id)):
                    print("-- warning : Particle ", str(id),"was not found at time: ", key)
                    continue

                body = self.timeLevels[key].getParticle(id)
                if(np.linalg.norm(body.force) == 0):
                    continue
                row_values =  ["{:<10}".format(str(key))] 
                row_values += ["{:<20.6f}".format(float(body.pos[i])) for i in range(3)]
                row_values += ["{:<20.6f}".format(float(body.linVel[i])) for i in range(3)]
                row_values += ["{:<20.6f}".format(float(body.force[i])) for i in range(3)]
                row_values += ["{:<20.6f}".format(float(np.linalg.norm(body.omega)))]
                row_values += ["{:<20.6f}".format(float(body.r))]
                row_values += ["{:<20.6f}".format(float(body.mass))]
                row_values += ["{:<20.6f}".format(float(body.kineticEnergy))]  
                file.write("\t".join(row_values) + "\n")        
    def calculateTotalKineticEnergy(self):
        with open("total_kinetic_energy.dat", "w") as file:
            file.write("{:<10}\t{:<20}\n".format("Time", "TotalKineticEnergy"))

            for time in self.getTimes():
                totalKineticEnergy = 0

                for particleID in self.timeLevels[time].getParticleList():
                    particle = self.timeLevels[time].getParticle(particleID)
                    particle.calculateKineticEnergy() 
                    totalKineticEnergy += particle.kineticEnergy 

                file.write("{:<10.6f}\t{:<20.6f}\n".format(time, totalKineticEnergy))

class timeLevel:
    def __init__(self):
        "Class to store the data of a sphere body."
        self.time = 0
        self.particles = {}

    def addParticle(self, particle):
        self.particles[particle.id] = particle
    
    def getParticleList(self):
        return list(self.particles.keys())

    def isParticlePresent(self, id):
        return id in self.particles

    def getParticle(self, id):
        return self.particles[id]

    def insertParticle(self, particle):
        self.particles[particle.id] = particle

    def exportParticlesToOF(self,outputDir):
        "Exports the particles to OpenFOAM format."
        particleList = self.getParticleList()
        particleList.sort(key = int)
        dirs = outputDir.split("/")
        outPutDir = ""
        for directory in dirs:
            outPutDir += directory+"/"
            if(not os.path.exists(outPutDir)):
                os.system("mkdir "+outPutDir)
        
        for id in particleList:
            self.particles[id].writeBodyInfoDict(outPutDir,"sphereFalling")

class bodyInfo:
    def __init__(self):
        "Class to store the data of a sphere body."
        self.id = 0
        self.pos = np.zeros(3)
        self.linVel = np.zeros(3)
        self.omega = np.zeros(3)
        self.r = 0
        self.mass = 0
        self.force = np.zeros(3)
        self.axis = np.ones(3)/np.linalg.norm(np.ones(3))
        self.kineticEnergy = 0


    def printStats(self):
        "Prints the body's stats. -- debuging purposes."
        print("ID: ",self.id)
        print("Position: ",self.pos)
        print("Linear Velocity: ",self.linVel)
        print("Angular Velocity: ",self.omega)
        print("Radius: ",self.r)
        print("Static: ",self.static)

    def calculateKineticEnergy(self):
        velocity_magnitude = np.linalg.norm(self.linVel)  # Compute magnitude of the velocity vector
        self.kineticEnergy = 0.5 * self.mass * velocity_magnitude**2

    def writeBodyInfoDict(self,path,particleName):
        "Writes the body info dictionary for OpenFOAM."
        if(np.linalg.norm(self.omega) != 0):
            self.axis =  self.omega/np.linalg.norm(self.omega)
        
        with open(str(path)+'body'+str(self.id)+'.info','w') as file:
            file.write("/*--------------------------------*- C++ -*----------------------------------*\\\n")
            file.write("  =========                 |\n")
            file.write("  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox\n")
            file.write("   \\\\    /   O peration     | Website:  https://openfoam.org\n")
            file.write("    \\\\  /    A nd           | Version:  8\n")
            file.write("     \\\\/     M anipulation  |\n")
            file.write("\\*---------------------------------------------------------------------------*/\n")
            file.write("FoamFile\n")
            file.write("{\n")
            file.write("    version     2.0;\n")
            file.write("    format      ascii;\n")
            file.write("    class       dictionary;\n")
            file.write("    location    \"/stonefly2/studenio/PhD_VirtMesh/FallingSpheres/bodiesInfo/0.005\";\n")
            file.write("    object      body"+str(self.id)+".info;\n")
            file.write("}\n")
            file.write("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n")
            file.write("\n")
            file.write("bodyId          "+str(self.id)+";\n")
            file.write("\n")
            file.write("bodyName        "+str(particleName)+";\n")
            file.write("\n")
            file.write("Vel             ( "+str(self.linVel[0])+" "+str(self.linVel[1])+" "+str(self.linVel[2])+" );\n")
            file.write("\n")
            file.write("omega           "+str(np.linalg.norm(self.omega))+";\n")
            file.write("\n")
            file.write("Axis            ( "+str(self.axis[0])+" "+str(self.axis[1])+" "+str(self.axis[2])+" );\n")
            file.write("\n")
            file.write("static          0;\n")
            file.write("\n")
            file.write("timeStepsInContWStatic 0;")
            file.write("\n")
            file.write("sphere\n")
            file.write("{\n")
            file.write("    position    ( "+str(self.pos[0])+" "+str(self.pos[1])+" "+str(self.pos[2])+" );\n")
            file.write("    radius      "+str(self.r)+";\n")
            file.write("}\n")
            file.write("\n")        
        file.close()

class bodeTimeResults:
    def __init__(self):
        "Class to store the data of a sphere body."
        self.id = 0
        self.time = []
        self.pos = []
        self.linVel = []
        self.omega = []
        self.force = []
        self.Impuls = []
        self.sumImpuls = 0
        self.traveledDistance = []
        self.sumDistance = 0
    
    def addTimeLevel(self,time, pos, linVel, omega, force):
        self.time.append(time)
        self.pos.append(pos)
        self.linVel.append(linVel)
        self.omega.append(omega)
        self.force.append(force)

    def calculateImpuls(self):
        self.sumImpuls = 0
        self.Impuls = [0]
        for i in range(0,len(self.force)-1):
            self.sumImpuls += (self.force[i]+self.force[i+1])*0.5*(self.time[i+1]-self.time[i])
            self.Impuls.append(I)
    
    def calculateTraveledDistance(self):
        self.sumDistance = 0
        self.traveledDistance = [0]
        for i in range(0,len(self.pos)-1):
            self.sumDistance += np.linalg.norm(self.pos[i+1]-self.pos[i])
            self.traveledDistance.append(self.sumDistance)


def readLiggghtsVTK(fileName,time):
    "Reads a VTK file and returns the data."
    LookedForStrings = ["POINTS ","VERTICES","id","type","mass ","v ","f ","radius ","omega "]
    Number_of_Points = []
    OuterList = [[],[],[],[],[],[],[]]
    ReadPairs = [[2,3],[4,5],[5,6],[6,7],[7,8],[8,6]]
    ContinueRequirement = False
    with open(fileName, 'r') as file:
        data = file.readlines()
    j = 0
    for i in range(len(data)):
        fInd = data[i].find(LookedForStrings[0]) #Look For Points
        if fInd == 0:
            Number_of_Points.append([float(strToVal) for strToVal in re.findall(r"-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?", data[i][fInd+len(LookedForStrings[0])::])])
            Number_of_Points_List  = [val[0] for val in Number_of_Points]
            Number_of_Points = Number_of_Points_List[0]
            if int(Number_of_Points) != 0:
                ContinueRequirement = True 

            if ContinueRequirement:
                for l in range(1,len(data)):
                    if data[i+l].find(LookedForStrings[1]) != 0 and ContinueRequirement: #Look For Vertricies
                        OuterList[0].append([float(strToVal) for strToVal in re.findall(r"-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?", data[i+l][::])])
                    elif data[i+l].find(LookedForStrings[1]) == 0 and ContinueRequirement: #Look For Vertricies
                        break            
        if ContinueRequirement:    
                readFirst = ReadPairs[j][0]
                readSecond = ReadPairs[j][1]
                fInd = data[i].find(LookedForStrings[readFirst])
                if fInd == 0 and ContinueRequirement:
                    k = 1
                    while k +i  < len(data):
                        if data[i+k].find(LookedForStrings[readSecond]) != 0:
                            OuterList[j+1].append([float(strToVal) for strToVal in re.findall(r"-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?", data[i+k][fInd::])])
                        elif data[i+k].find(LookedForStrings[readSecond]) == 0 or i+k == len(data):
                            j+=1
                            break
                        k += 1
    #Read the data
    #Process the data
    bodyPos = []
    ParticleIDs = []
    Mass_Particles = []
    v_Particles = []
    f_Particles = []
    radius_Particles = []
    omega_Particles = []
    TimeDependentLists = [bodyPos,ParticleIDs,Mass_Particles,v_Particles,f_Particles,radius_Particles,omega_Particles]
    #Merges All SubList into one
    for i in range(0,len(OuterList)):
        for j in range(0,len(OuterList[i])):
            for k in range(0,len(OuterList[i][j])):
                val = (OuterList[i][j][k])
                TimeDependentLists[i].append(val)
    #Merges All SubList into one
    bodyPos = TimeDependentLists[0]
    ParticleIDs = TimeDependentLists[1]
    Mass_Particles = TimeDependentLists[2]
    v_Particles = TimeDependentLists[3]
    f_Particles = TimeDependentLists[4]
    radius_Particles = TimeDependentLists[5]
    omega_Particles = TimeDependentLists[6]

    currTimeLevel = timeLevel()
    timeLevel.time = time

    for i in range(0,len(ParticleIDs)):
        body = bodyInfo()
        body.id = int(ParticleIDs[i])
        body.pos = np.array([bodyPos[3*i],bodyPos[3*i+1],bodyPos[3*i+2]])
        body.linVel = np.array([v_Particles[3*i],v_Particles[3*i+1],v_Particles[3*i+2]])
        body.omega = np.array([omega_Particles[3*i],omega_Particles[3*i+1],omega_Particles[3*i+2]])
        body.r = radius_Particles[i]
        body.mass = Mass_Particles[i]
        body.force = np.array([f_Particles[3*i],f_Particles[3*i+1],f_Particles[3*i+2]])

        currTimeLevel.insertParticle(body)

    return currTimeLevel

def main():
    "Main function"
    #==========================#
    dt = 1e-6
    nDigits = 2
    inputName = "particles_Res_"
    dirName = "post/"
    #==========================#
    presentResults = getTimeLevelList(dirName,inputName)
    # presentResults= [presentResults[0]]
    simRes = simulationResults()
    for t in presentResults:
        time = t*dt
        time = round(time,nDigits)
        fileName = dirName+inputName+str(t)+".vtk"
        timeLevel = readLiggghtsVTK(fileName,time)
        simRes.addTimeLevel(timeLevel,time)
        simRes.addBodiesToGlobalList()
        print("simRes.globalIDList: ", simRes.globalIDList)

    simRes.calculateTotalKineticEnergy()
    
    #Convert initial stateTo OpenFOAM format
    tList = simRes.getTimes()
    tList.sort(key = float)
    i = 0

    for id in simRes.globalIDList:
        simRes.getBodyResults(id)
        simRes.getBodyNonZeroForces(id)
    
if __name__ == "__main__":
    main()