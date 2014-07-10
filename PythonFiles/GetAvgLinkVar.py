from AAPI import *
from angParams import ResultFileName
from secIDs_table import secIDs
from secSourceIDs_table import secSourceIDs
from AllSections import AllSections
import sys
import csv
# the file angParams.py is in the TR/PyInputs directory
# we only record events after the warmupPeriod (900sec ie 15min)
matrixFilePath = 'C:/Users/Krishna/Dropbox/CriticalOD/PythonFiles/ODpairsListForResult.txt'
#global variables
file2 = 0
totNbVeh = 0
totTT = 0 

def AAPILoad():
	return 0

def AAPIInit():
	#AKIPrintString( "AAPIInit" )
	return 0

def AAPIManage(time, timeSta, timeTrans, acycle):
	return 0

def AAPIPostManage(time, timeSta, timeTrans, acycle):
	return 0

def AAPIFinish():
	global file2
	global totExtraVeh
	#AKIPrintString( "AAPIFinish" )

	#Get the average link travel time deviation
	nbSecs = len(AllSections)
	SumTravelTimeDeviation = 0
	SumLinkVehicles = 0
	for secIter in range(0,nbSecs):
		currSecID = AllSections[secIter] 
		sectStruct = AKIEstGetGlobalStatisticsSection(currSecID,0)
		currTime = sectStruct.TTd
		currNbVeh = sectStruct.Flow
		SumTravelTimeDeviation = SumTravelTimeDeviation + currTime*currNbVeh
		SumLinkVehicles = SumLinkVehicles + currNbVeh;
	AvgLinkTTDeviation = (SumTravelTimeDeviation/60)/(SumLinkVehicles)

	crs = open(matrixFilePath, "r")
	iter = 0;
	for columns in ( raw.strip().split() for raw in crs ):  
		iter = iter+1;
	matrix = [[0]*2 for i in range(iter)]
	crs.closed
	crs = open(matrixFilePath, "r")

	SumODTravelTimeDeviation = 0
	SumODLinkVehicles = 0
	iter = 0;
	for columns in ( raw.strip().split() for raw in crs ):
		matrix[iter][0]=int(columns[0])
		matrix[iter][1]=int(columns[1])
		iter = iter+1;
	crs.closed

	#for ODIter in range(0,iter):
	ODStruct = AKIEstGetGlobalStatisticsODPair(3695,3753, 0)
	#SumODTravelTimeDeviation = SumODTravelTimeDeviation + ODStruct.TTa
	#SumODLinkVehicles = SumODLinkVehicles + ODStruct.Flow
	#AvgODTTDeviation = (SumODTravelTimeDeviation/60)/(iter)

	file2 = open(ResultFileName, 'a')
	file2.write('%f %f %f %f %d\n'%(totTT/(totNbVeh*60),0,AvgLinkTTDeviation,ODStruct.TTa,ODStruct.report))
	file2.close()

	return 0

def AAPIUnLoad():
	return 0
	
def AAPIPreRouteChoiceCalculation(time, timeSta):
	return 0


############################################################################

def AAPIEnterVehicle( a_nVehId, a_nSectionId ):
	AKIVehSetAsTracked( a_nVehId )
	return 0


############################################################################

def AAPIExitVehicle( a_nVehId, a_nSectionId ):
	global totNbVeh
	global totTT
	if (AKIGetCurrentSimulationTime()>900):
		infVeh = AKIVehTrackedGetInf( a_nVehId )
		AKIVehSetAsNoTracked( a_nVehId )
		totNbVeh = totNbVeh+1
		totTT = totTT + (AKIGetCurrentSimulationTime()-infVeh.SystemEntranceT)
	return 0



############################################################################
####

