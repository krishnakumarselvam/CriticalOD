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


#This script creates a new aimsun file, with the new OD demand matrix = old +1
#It also writes the old OD matrix into a text file

from PyANGBasic import *
from PyANGKernel import *
from PyANGConsole import *
from PyANGAimsun import *
import sys

from angParams import newAngName
from angParams import baseAngName
from angParams import textFileName
matrixFilePath = 'C:/Users/Krishna/Dropbox/CriticalOD/Aimsun/ODpairs.txt'
def main( argv=None):
        
        if argv is None:
           argv = sys.argv
        
        console = ANGConsole()
        print "Program Compiled"
        initModel = baseAngName
       
        if console.open( initModel ):
                file = open( matrixFilePath, 'w' )
                objectType = GKSystem.getSystem().getActiveModel().getType( "GKODMatrix" )
                #centroids = console.getModel().getCatalog().getObjectsByType( objectType )
                for currODMatrix in console.getModel().getCatalog().getObjectsByType( objectType ).itervalues():
                        temp = currODMatrix.getId()
                        if temp==3377 :
                                print temp
                                centroids = currODMatrix.getCentroidConfiguration().getCentroidsInOrder()
                                #now we have a list of all the centroids
                                for origin in centroids:
                                        for destination in centroids:
                                                if origin != destination:
                                                        trips = currODMatrix.getTrips( origin, destination)
                                                        if trips > 0:
                                                                file.write( '%u %u %f\n' % (origin.getId(), destination.getId(), trips)  )
                                                                currODMatrix.setTrips( origin, destination, trips+1)
                               
     
                newReplic = GKSystem.getSystem().newObject("GKReplication", console.getModel())
                currExp = console.getModel().getCatalog().find( int( 3370 ) )
                #currSeed = random.randrange(0,500000,1)
                #newReplic.setRandomSeed(currSeed)
                currExp.addReplication(newReplic)
                print "*******************************************"
                print "Model Created with no errors"
                print "*******************************************"
        console.save( newAngName )     
        console.close()
        GKSystem.deleteSystem()
        
if __name__ == "__main__":
        sys.exit( main() )
