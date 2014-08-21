from AAPI import *
from angParams import newAvgTTFile
from secIDs_table import secIDs
from secSourceIDs_table import secSourceIDs
from tradeoffSD import tradeoff
# the file angParams.py is in the TR/PyInputs directory
# we only record events after the warmupPeriod (900sec ie 15min)

#global variables
file2 = 0
totNbVeh = 0
totTT = 0 
totExtraVeh = 0

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

	nbSecs = len(secIDs)
	totTTAllSect = 0
	totVarAllSect = 0
	for secIter in range(0,nbSecs):
		currSecID = secIDs[secIter] 
		sectStruct = AKIEstGetGlobalStatisticsSection(currSecID,0)
		currTime = sectStruct.TTa
		currSD = sectStruct.TTd
		currVar = currSD**2
		# total travel time is not cannot be accessed, so I use TTa which is a rounded version of TotalTravelTime/Flow for each section
		# displayed in seconds
		currNbVeh = sectStruct.Flow
		# multiplied by flow so its the total time (not averaged link time)
		# so that final average weights different links according to their flow
		totTTAllSect = totTTAllSect + currTime/60
		totVarAllSect = totVarAllSect + (currVar/3600)
		

	#nbSecsSource = len(secSourceIDs)
	#totEntryVeh = 0
	#for secIter in range(0,nbSecsSource):
	#	currSecID = secSourceIDs[secIter] 
	#	sectStruct = AKIEstGetGlobalStatisticsSection(currSecID,0)
	#	totEntryVeh = totEntryVeh + sectStruct.Flow

	# [min]
	#avgSubNwkTT = (totTTAllSect/60)/(totEntryVeh+totExtraVeh)
	totSDAllSect = totVarAllSect**0.5
	TTSD = totTTAllSect+tradeoff*totSDAllSect
	file2 = open(newAvgTTFile, 'a')
	# write  totExtrVeh[min], totNbVeh, avgTT[min],  replicID
	file2.write('%f %f %f %d\n'%(totTTAllSect,totSDAllSect,TTSD,ANGConnGetReplicationId()))
	file2.close()

	return 0

def AAPIUnLoad():
	return 0
	
def AAPIPreRouteChoiceCalculation(time, timeSta):
	return 0


############################################################################

def AAPIEnterVehicle( a_nVehId, a_nSectionId ):
 #   global totExtraVeh

 #   if (a_nSectionId==142):
 #      totExtraVeh = totExtraVeh+1
 #   if (a_nSectionId==759):
 #      totExtraVeh = totExtraVeh+1
 #   if (a_nSectionId==760):
 #      totExtraVeh = totExtraVeh+1

    return 0


############################################################################

def AAPIExitVehicle( a_nVehId, a_nSectionId ):
    return 0


############################################################################
####

