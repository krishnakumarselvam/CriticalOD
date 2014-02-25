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


#This script creates a new aimsun file, with the new OD demand matrix

from PyANGBasic import *
from PyANGKernel import *
from PyANGConsole import *
from PyANGAimsun import *
import sys
import csv
from angParams import newAngName
from angParams import baseAngName
from angParams import ResultFileName
from angParams import ODFileName
matrixFilePath = 'C:/Users/Krishna/Dropbox/CriticalOD/Aimsun/ODpairs.txt'
def main( argv=None):
        
        if argv is None:
           argv = sys.argv
        
        console = ANGConsole()
        #print "Program Compiled"
        initModel = baseAngName
        #print ODFileName
       
        if console.open( initModel ):

                #Read inputs
                crs = open(ODFileName, "r")
                iter = 0;
                for columns in ( raw.strip().split() for raw in crs ):  
                    iter = iter+1;
                matrix = [[0]*3 for i in range(iter)]
                crs.closed
                
                crs = open(ODFileName, "r")
                iter = 0;
                for columns in ( raw.strip().split() for raw in crs ):
                    matrix[iter][0]=int(columns[0])
                    matrix[iter][1]=int(columns[1])
                    matrix[iter][2]=float(columns[2])
                    iter = iter+1;
                crs.closed
                ctr=0

                file = open( matrixFilePath, 'w' )
                objectType = GKSystem.getSystem().getActiveModel().getType( "GKODMatrix" )
                #centroids = console.getModel().getCatalog().getObjectsByType( objectType )
                for currODMatrix in console.getModel().getCatalog().getObjectsByType( objectType ).itervalues():
                        temp = currODMatrix.getId()
                        if temp==3377 :
                                #print temp
                                centroids = currODMatrix.getCentroidConfiguration().getCentroidsInOrder()
                                #now we have a list of all the centroids
                                for origin in centroids:
                                        for destination in centroids:
                                                if origin != destination:
                                                        trips = currODMatrix.getTrips( origin, destination)
                                                        if trips > 0:
                                                                currODMatrix.setTrips( origin, destination, matrix[ctr][2])
                                                                ctr=ctr+1;
                                                                                                                        
                               
                               

                newReplic = GKSystem.getSystem().newObject("GKReplication", console.getModel())
                currExp = console.getModel().getCatalog().find( int( 3370 ) )
                #currSeed = random.randrange(0,500000,1)
                #newReplic.setRandomSeed(currSeed)
                currExp.addReplication(newReplic)
                #print "*******************************************"
                #print "Model Created with no errors"
                #print "*******************************************"
                print"..."
                              
        console.save( newAngName )     
        console.close()
        GKSystem.deleteSystem()
 
        
if __name__ == "__main__":
        sys.exit( main() )
