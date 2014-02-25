# differs from  simRepAll_newMdl.py : the seed used is currReplicID 
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
from angParams import newAngName
# angParams.py are in the TR/PyInputs directory

# Simulates a replication. It ask to thr AIMSUN Micro plugin (internal name:
# GGetram) for a simulator, the adds a simuation task and simulates it
def simulate( replication ):
	plugin = GKSystem.getSystem().getPlugin( "GGetram" )
	if plugin != None:
		simulator = plugin.getCreateSimulator( replication.getModel() )
		if simulator != None:
			simulator.addSimulationTask( replication, GKReplication.eBatch)
			simulator.simulate()

# Loads a network (newAngName) and simulates nbReplic replications on it
#
def main( argv=None):

	if argv is None:
		argv = sys.argv
	
	console = ANGConsole()
	print newAngName
	if console.open( newAngName ):
	   	print 'model to be run:'
		print newAngName
		replicType = GKSystem.getSystem().getActiveModel().getType( "GKReplication" )
		# check that there is only one replication defined
		j = 0
		for currRep in console.getModel().getCatalog().getObjectsByType( replicType ).itervalues():
			j=j+1
		if (j==1):
			for currRep in console.getModel().getCatalog().getObjectsByType( replicType ).itervalues():
				if (currRep.getSimulationStatus()>1):
					print 'replic is neither pending nor done. check status.'

				if currRep != None:
					simulate( currRep )
					print 'replication done:'
					print currRep.getId()
		else:
			print 'more than 1 replic in model, so we didnt launch aimsun'
	console.close()
	GKSystem.deleteSystem()
		
# Entry code, the script starts here:
if __name__ == "__main__":
	sys.exit( main() )
