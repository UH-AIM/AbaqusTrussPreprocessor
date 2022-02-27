#load modules
from abaqus import *
from abaqusConstants import *
backwardCompatibility.setValues(includeDeprecated=True, reportDeprecated=False)
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
import mesh
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
import regionToolset
import getpass

import csv
import numpy as np
import math
import os
import sys

src_file_path = inspect.getfile(lambda: None)
fDir = src_file_path.rsplit('/', 1)[0]

# The main function should contain all user alterable inputs
def main():
	# fPath is where I store all of my benchmark models
	fPath=fDir+'/input'
	# This is the subfolder with the name of the model
	fName='25_bar'
	# This is where you want where all the analysis files to go
	wDir=fDir+'/output/'+fName
	# If the folder doesn't exist, it will be created
	if not os.path.exists(wDir):
		os.makedirs(wDir)
	# because of scaling issues, sometimes I use mm, but mostly m
	unit='m'
	# there could be multiple load cases
	case=1
	session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)
	# Constructor reads in all the files
	objAbaModel=cAbaModel(wDir,fPath,fName,unit,case)
	# This is where the geometry is generated
	objAbaModel.importGeom()
	# defAbaMat function has the constitutive model of all connex materials
	objAbaModel.defAbaMat('elastic')
	# Read cross section from file
	objAbaModel.crossSection()
	# Create load step
	objAbaModel.loadStep()
	# Assignment of boundary conditions
	objAbaModel.BC()
	# create a job
	objAbaModel.abaJob()
	# run the job
	# objAbaModel.abaJobSubmit()
	
class cAbaModel(object):
	def __init__(self, wDir,fPath,fName,unit,case):
		os.chdir(wDir)
		
		listMatConnex=['polymer']
		Mdb()
		
		mdb.models.changeKey(fromName='Model-1', toName=fName)
		self.m=mdb.models[fName]
		self.dUnit=1 # if mm, then this is 1000 i think
		self.fName=fName
		self.fInName=fPath+'/'+fName+'/'+unit+'/'
		self.partName='part-'+fName
		self.stepName='step-'+fName
		self.instanceName='instance-'+fName
		self.jobName='job-'+fName
		
		self.case=case
		# Element connectivity file is called con_#.csv
		reader=list(csv.reader(open(self.fInName+'con_'+str(self.case)+'.csv','rb'),delimiter=','))
		# read in node start and node end
		self.listElem=[j[0:2] for j in reader]
		self.listElem=np.array(self.listElem).astype('int')
		# read in cross section area (for truss, shape doesn't matter)
		self.listArea=[j[2] for j in reader]
		self.listArea=np.array(self.listArea).astype('float')
		self.listMat=[j[3] for j in reader]
		# read in material index (1 - 14)
		self.listMat=np.array(self.listMat).astype('str')
		# Read in coordinate list and put into listNodes
		reader=list(csv.reader(open(self.fInName+'coord_'+str(self.case)+'.csv','rb'),delimiter=','))
		self.listNodes=[r[0:3] for r in reader]
		self.listNodes=np.array(self.listNodes).astype('float')
		# read in all the forces
		reader=list(csv.reader(open(self.fInName+'forces_'+str(self.case)+'.csv','rb'),delimiter=','))
		self.listForces=[r[0:3] for r in reader]
		self.listForces=np.array(self.listForces).astype('float')
		# read in all the boundary conditions
		reader=list(csv.reader(open(self.fInName+'constraints_'+str(self.case)+'.csv','rb'),delimiter=','))
		self.listConstraints=[r[0:3] for r in reader]
		self.listConstraints=np.array(self.listConstraints).astype('float')
		
	def importGeom(self):
		# create nodes 
		p=self.m.Part(self.partName)
		listObjNode=[];
		self.listObjElem=[];
		for i,coordNodes in enumerate(self.listNodes):
			objNode=p.Node(coordNodes)
			listObjNode.append(objNode)
		# create elements
		for j,pairElem in enumerate(self.listElem):
			objElem=p.Element([listObjNode[pairElem[0]-1],listObjNode[pairElem[1]-1]],LINE2)
			self.listObjElem.append(objElem)
		
	def crossSection(self):
		p = self.m.parts[self.partName]
		e = p.elements
		for i in range(0,len(self.listElem)):
			elements = e[i:i+1]
			region = regionToolset.Region(elements=elements)
			
			self.m.TrussSection(name='Sec-'+str(i), material=self.listMat[i], area=self.listArea[i])
			p.SectionAssignment(region=region, sectionName='Sec-'+str(i), offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)
		elements = e[:]
		region = regionToolset.Region(elements=elements)
		orientation=None
		p.MaterialOrientation(region=region, orientationType=GLOBAL, axis=AXIS_1, additionalRotationType=ROTATION_NONE, localCsys=None, fieldName='', stackDirection=STACK_3)
		
		elemTrussType = mesh.ElemType(elemCode=T3D2, elemLibrary=STANDARD)
		
		pickedRegions =(elements, )
		p.setElementType(regions=pickedRegions, elemTypes=(elemTrussType, ))
		p.regenerate()
		
		a = self.m.rootAssembly
		a.Instance(name=self.instanceName, part=p, dependent=ON)
		a.regenerate()
	def BC(self):
		a = self.m.rootAssembly
		n = a.instances[self.instanceName].nodes
		
		for i,x in enumerate(self.listConstraints):
			region = regionToolset.Region(nodes=n[int(x[0])-1:int(x[0])])
			if x[1]==1:
				bcList=[UNSET,UNSET,x[2]]
			elif x[1]==2:
				bcList=[UNSET,x[2],UNSET]
			elif x[1]==4:
				bcList=[x[2],UNSET,UNSET]
				
			self.m.DisplacementBC(name='BC-'+str(i), 
			createStepName=self.stepName, region=region, u1=bcList[0], u2=bcList[1], u3=bcList[2], ur1=UNSET, ur2=UNSET, ur3=UNSET, amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)
			
		for i,x in enumerate(self.listForces):
			region = regionToolset.Region(nodes=n[int(x[0])-1:int(x[0])])
			if int(x[1])==1:
				lcList=[UNSET,UNSET,x[2]]
			elif int(x[1])==2:
				lcList=[UNSET,x[2],UNSET]
			elif int(x[1])==4:
				lcList=[x[2],UNSET,UNSET]
				
			self.m.ConcentratedForce(name='Load-'+str(i), createStepName=self.stepName, region=region, cf1=lcList[0], cf2=lcList[1], cf3=lcList[2], distributionType=UNIFORM, field='', localCsys=None)

	def loadStep(self):
		self.m.StaticStep(name=self.stepName, previous='Initial', initialInc=0.01, minInc=1e-30, maxInc=0.05, nlgeom=ON)
		
	def abaJob(self):
		mdb.Job(name=self.jobName, model=self.fName, description='', type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=1, numGPUs=0)
		
	def abaJobSubmit(self):
		mdb.jobs[self.jobName].submit(consistencyChecking=OFF)
		mdb.jobs[self.jobName].waitForCompletion()
    
	def defAbaMat(self, matModel):
		listMat=[
			{'name':'polymer','density':1180,'E':2558647069.0, 'nu':0.4}]

		for i, dicMat in enumerate(listMat):
			self.m.Material(name=dicMat['name'])
			self.m.materials[dicMat['name']].Density(table=((dicMat['density']*self.dUnit**4,),))
			if matModel=='elastic':
				self.m.materials[dicMat['name']].Elastic(table=((dicMat['E']*self.dUnit**2, dicMat['nu']), ))
                
	def fMC(self,cCoord1,cCoord2): # findMeanCoord
		meanCoord=np.zeros([1,3])
		meanCoord[0,0:2]=(np.array(cCoord1) - np.array(cCoord2))/2 + np.array(cCoord2)
		meanCoord=meanCoord[0]
		return meanCoord
		
	def expC(self,cCoord):
		meanCoord=np.zeros([1,3])
		meanCoord[0,0:2]=cCoord
		meanCoord=meanCoord[0]
		return meanCoord

if __name__=='__main__':
	main()
