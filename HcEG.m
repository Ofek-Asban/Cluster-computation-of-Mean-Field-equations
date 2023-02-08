function [] = SlaveHC2(RealNum, N) % This SlaveHeatingAndCooling file gets a number of realization, lowerst and highest temps and then heating and cooling the sample between thous temps in 20 steps and then saving the E vector into a file with the specific realization number RealNum - number of specific realization

rng('shuffle') %This command in each matlab run change the seed to produce different random numbers (Only relevant on distributed computing because each jobs that goes to new porcessor opens matlab for the first time which gives same seed that generetes same random numbers)

T_low = 0.1; 

%%%Initial temperature which will be the final temperature at the end of the process (0.05 in Ariel's paper for calculation of the DOS and 0.1 for calculation of the distribution of eigen values and in case of weak interaction e^2/(T*r_nn)=1 and mot 20 T will be 1 and therefore W=10 the disorder energy where W/T=10)

T_high = 10; 

%%%hottest temperature (in the case of strong interactions T_high=1 bcs there the temp will be comperable with the interactions though e^2/(T*r_nn)=20)

W=1; %disorder energy

epsilon = W*(rand(1,N)-0.5); %calibrating random on-site energies in the range [-W/2, W/2] where W is th disorder energy which is W/T=10

E = 20*(rand(1,N)-0.5); % calibrating empty vertor for input the average site energies

[R_ij, r_nn] = DistanceMatrix(N); %creating new distance matrix in each realization because this gives the same DOS in low number of sites (it will help to get real uniform distribution of distances) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%Heating and then Cooling procedure%%%%%%%%%%%%%%%%%%%%%%%%%% 
T = T_low; % initializing the temperature
Num_Of_T_Steps = 10; %number of heating and then cooling steps (ref 18 used 20 cooling steps)
size_Of_T_Step = 2*(T_high-T_low)/(Num_Of_T_Steps); %size of each step (multiplyed by 2 because we need to heat and then cool down the system)

for k = 1:(Num_Of_T_Steps+1)
    
    [E] = EnergyIterationsSlave2(N, epsilon, E, R_ij, T, r_nn); %calling the function to get the solution of E vector for a given T and saving it in a file
        
    if(k <= Num_Of_T_Steps/2) %Heating the system
        T = T + size_Of_T_Step;
    end
    
    if(k > Num_Of_T_Steps/2) %Cooling back the system
        T = T - size_Of_T_Step;
    end	
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd /fastspace/users/asban/output_mat %changing again the directory to output_mat to sae there the SleaveHeatingAndCooling_RealNum file
fname = sprintf('DmatrixAndE_%d.mat', RealNum); %adding the realization number to the file name
save(fname, 'R_ij', 'r_nn', 'E'); %saving the E vector to a new file of current Realization (these E vector in the file will represent a realization of the system after heating and cooling the system).

