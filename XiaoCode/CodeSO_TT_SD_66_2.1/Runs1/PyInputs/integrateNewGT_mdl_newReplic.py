################################################################################
#                                                                              #
# (c) 2006 TSS-Transport Simulation Systems                                    #
#                                                                              #
################################################################################

# This is a script for the AIMSUN Console application. It loads a network (from
# the first argument) and simulates a replication on it (the second argument of
# the script).
#
# The execution line is:
# angconsole.exe -script simulate.py my_network.ang 100
# Where:
# - my_network.ang is the name of the network file to load. It can contain the
#   full network path as "c:/Networks/my_network.ang"
# - 100 is the replication id, change it with a valid one

# Disclaimer of Warranty: TO THE EXTENT PERMITTED BY APPLICABLE LAW THE AIMSUN
# SCRIPTS, INCLUDING BUT NOT LIMITED TO ALL DATA, TOOLS AND CALCULATIONS THEREIN
# ARE PROVIDED "AS IS" AND WITHOUT EXPRESS OR IMPLIED WARRANTY OF ANY KIND BY
# EITHER TSS OR ANYONE ELSE WHO HAS BEEN INVOLVED IN THE CREATION, PRODUCTION OR
# DELIVERY OF THE AIMSUN SCRIPTS, INCLUDING BUT NOT LIMITED TO ANY EXPRESS OR
# IMPLIED WARRANTY OF MERCHANTABILITY, ACCURACY, QUIET ENJOYMENT,
# NONINFRINGEMENT OR FITNESS FOR A PARTICULAR PURPOSE. NO COVENANTS, WARRANTIES
# OR INDEMNITIES OF ANY KIND ARE GRANTED BY TSS TO YOU THE USER. SHOULD THE
# PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING,
# REPAIR OR CORRECTION.

from PyANGBasic import *
from PyANGKernel import *
from PyANGConsole import *
from PyANGAimsun import *
import sys
from generic_new_gtContDico import controlDico
from angParams import newAngName
# generic_new_gtContDico  and angParams.py are in the TR2/PyInputs directory


def main( argv=None):

	if argv is None:
	   argv = sys.argv
	
	console = ANGConsole()
	initModel = 'C:/MIT/Lab/KanchanaNanduri/CodeSO_TT_SD/Runs1/Aimsun/CentreV5_getAvgTT_NoReplic.ang'
    # N:/RunsTR/Quad30/Runs15_subnwkTT/Aimsun/CentreV5_getAvgTT_NoReplic.ang'
	if console.open( initModel ):
		#if console.open( argv[1] ):
	   #fileNewOutMdl = 'C:/Program Files/TSS-Transport Simulation Systems/AIMSUN_simLaus_docs/tmpNewGeneratedGT.ang'
	   
	   #print argv[2]
	   #from argv[2] import controlDico 
	   print controlDico

	   controlPlanType = GKSystem.getSystem().getActiveModel().getType( "GKControlPlan" )
	   #controlPlanType = GKSystem.getSystem().getType( "GKControlPlan" )
	   for currControlPlan in console.getModel().getCatalog().getObjectsByType( controlPlanType ).itervalues():

	       fileName = 'C:/MIT/Lab/KanchanaNanduri/CodeSO_TT_SD/Runs1/GTsUsed/diffPlans_aimsunOutput.txt'
	       # fileName = 'N:/RunsTR/Quad30/Runs15_subnwkTT/GTsUsed/diffPlans_aimsunOutput.txt'
	       file = open(fileName, 'w')
	       file.write('Old vs. new control plan signals\n')
	       file.write('OLD: nodeID; phaseID; phaseDuration; controlPlanSignalId; turningID; originSection; fromLane; toLane; destSection; greenTimeDuration[sec];  totalCycleDuration[sec]; start[sec]; \n')		
	       file.write('NEW: nodeID; phaseID; phaseDuration; controlPlanSignalId; turningID; originSection; fromLane; toLane; destSection; greenTimeDuration[sec];  totalCycleDuration[sec]; start[sec]; \n\n')
	       	       
	       allControlJunctionsAndIds = currControlPlan.getControlJunctions()    
	       for currId, currJctn in allControlJunctionsAndIds.items():	
				if controlDico.has_key(currId):
				# # if (currId == nodeTest):
					currDico_phases = controlDico[currId]
					currControlPlanSignals = currJctn.getSignalsStartDuration() 
					allPhases = currJctn.getPhases()
					if (len(currControlPlanSignals)>0):
						# its a signalized intersection
						# for a given phase find the turnings: sections and lanes, that belong to it
						for currPhase in allPhases:
							if currDico_phases.has_key(currPhase.getId()):
								# # if (currPhase.getId() == phaseTest):
								currDico_thisPhase = currDico_phases[currPhase.getId()]
								# currDico_thisPhase contains: the difference in start times and the new green times of the current phase
								iterK = [1,2]
								for k in iterK:
									if (k ==1): 	
										file.write('Old:\n ')
									else: 
										file.write('New: \n')
										# change the start time if needed
										if (currDico_thisPhase[0] != 0):
											oldStart = currPhase.getFrom()
											newStart = oldStart + currDico_thisPhase[0] 
											#if (currId == 944):
											#print oldStart
											#print 'new'
											#print newStart
											# python indexes arrays starting at 0
											currPhase.setFrom(newStart)
					 
										# change the green time if needed
										oldGreen = currPhase.getDuration()
										newGreen = currDico_thisPhase[1] 
										if (oldGreen != newGreen):
											currPhase.setDuration(newGreen)

									phaseSignals = currPhase.getSignals()
									for currSignal in phaseSignals:
										# find the corresponding controlPlanSignal (currSignal is a controlPhaseSignal)
										for currControlPlan in currControlPlanSignals:
											currSig2 = currControlPlan.getControlPlanSignal()
											if (currSig2.getId() == currSignal.signal):
												thisTurns = currSig2.getTurnings()
												for  currTurn in thisTurns:
													file.write('%i'%(currId))
													file.write('	%i'%(currPhase.getId()))
													file.write(' 	%i'%(currPhase.getDuration()))				                           
													file.write('	%i'%(currSig2.getId()))
													file.write(' 	%i'%(currTurn.getId()))
													file.write(' 	%i'%(currTurn.getOrigin().getId()))
													file.write('	%i'%(currTurn.getOriginFromLane()))
													file.write('	%i'%(currTurn.getOriginToLane()))
													file.write(' 	%i'%(currTurn.getDestination().getId()))
													file.write(' 	%i'%(currControlPlan.getDuration()))		
													file.write(' 	%i'%(currJctn.getCycle()))
													file.write(' 	%i'%(currControlPlan.getStart()))		
													file.write('\n')	
					
	       file.close()


	   # create a new replication; add it to expt 3370; write the replication ID into the file
	   newReplic = GKSystem.getSystem().newObject("GKReplication", console.getModel())
	   currExp = console.getModel().getCatalog().find( int( 3370 ) )
	   #currSeed = random.randrange(0,500000,1)
	   #newReplic.setRandomSeed(currSeed)
	   currExp.addReplication(newReplic)


	   fileName2 = 'C:/MIT/Lab/KanchanaNanduri/CodeSO_TT_SD/Runs1/PyOutputs/currRepID.txt'
	   # fileName2 = 'N:/RunsTR/Quad30/Runs15_subnwkTT/PyOutputs/currRepID.txt'
	   file = open(fileName2, 'w')
	   file.write('%i'%(newReplic.getId()))
	   file.close()

	   fileName3 = 'C:/MIT/Lab/KanchanaNanduri/CodeSO_TT_SD/Runs1/PyOutputs/allMdls_repIDs_seeds.txt'
	   # fileName3 = 'N:/RunsTR/Quad30/Runs15_subnwkTT/PyOutputs/allMdls_repIDs_seeds.txt'
	   file = open(fileName3, 'a')
	   file.write('%i	'%(newReplic.getId()))
	   file.write('%i	'%(newReplic.getRandomSeed()))
	   file.write('%s	\n'%(newAngName))
	   file.close()


	console.save( newAngName )     
	console.close()
	
# Entry code, the script starts here:
if __name__ == "__main__":
	sys.exit( main() )
