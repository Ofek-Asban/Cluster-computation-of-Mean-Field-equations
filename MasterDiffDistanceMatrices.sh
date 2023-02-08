#!/bin/csh

@ narg = $#argv

# Checking with the user if the number and order of the input parameters is correct
if($narg != 3) then
        echo "usage: Master.sh <queue name> <Number Of Realizations> <N-number of sites>"
        exit
endif

# Setting parameters
set q = $1.q              #Getting the queue name
@ rNum = $2          #Getting the number of realizations 
@ N = $3                 #Number of sites
@ count = 1            #Initializing realization counter
@ countJOBS = 1  #Initializing the job counter

# Calling the HeatingAndCooling.m file rNum times. This makes rNum solutions of random realizations of energy E vectors (Each solution is obtrained after 20 steps of heating and cooling the system to optimized the solution) 
while ($count <= $rNum)

           #Calling HC.m and sending it 4 parameters and in turn the srcipt runs the SlaveHC.m file with these paramaters and creats a file SlaveHeatingAndCooling_$count.mat
           #count is the current realizations number. 
           qsub -cwd -q $q << EOJ
           cd /fastspace/users/asban/parallelTLS/storage/matlab/bin/matlab << M_PROG
           HcTLS(${count}, ${N});
           M_PROG 
           EOJ

           set countJOBS = `qstat | grep -c asban` # saving the number of jobs in the queue 
           @ countJOBS = $countJOBS # changing the number of JOBS to an interger
	
           #Limiting the number of running jobs to 1000 to aviod over use of resources.
           while ($countJOBS >= 1000)
	echo "$countJOBS in queue, waiting for free cores to continue with DmatrixAndE_$count+1"	
	sleep 30
	set countJOBS = `qstat | grep -c asban`  #Checking the current numer of jobs
	@ countJOBS = $countJOBS 	
            end
	
            @ count++
end

echo "Finished"
