# This script checks if some of the JOBS crashed during computation and run again these jobs. 
#!/bin/csh

set Path = "/fastspace/users/asban/output_mat"
set q = $1.q       #getting queue nume
@ rNum = $2   #getting realization number
@ N = $3          #number of sites
@ count = 1 

while($count <= $rNum)
        
       if (! -f $Path/Edelta0Dmat_$count.mat) then 

               qsub -cwd -q $q << EOJ
               cd /fastspace/users/asban/parallelTLS/storage/matlab/bin/matlab << M_PROG
               HcTLS(${count}, ${N});
               M_PROG
               EOJ

               echo "File Eanddelta0_$count.mat was crushed or override during computation, resending it again to queue."

               sleep 1 #keeping time between each queue so there will be no overload on the license server.

       endif

       @ count++
end

echo "finished going on $rNum files"
