from AAPI import *
from angParams import newAvgTTFile
# the file angParams.py is in the TR/PyInputs directory
# we only record events after the warmupPeriod (900sec ie 15min)

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
	global totNbVeh
	global totTT
	#AKIPrintString( "AAPIFinish" )


	# [min]
	file2 = open(newAvgTTFile, 'a')
    # write totTT [min], totNbVeh, avgTT[min], replicID
	file2.write('%f %d %f %d\n'%(totTT/60,totNbVeh,totTT/(totNbVeh*60),ANGConnGetReplicationId()))
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

