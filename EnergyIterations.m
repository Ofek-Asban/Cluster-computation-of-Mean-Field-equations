function [E, count] = EnergyIterationsSlave2(N, epsilon, E, R_ij, T, r_nn) %this function number of sites N calling the function distanceMatrixandRnn and calculated metastable solution to eq(5) and returns the solution in a shape of a vectors E for average site enegies

rng('shuffle') %this command in each matlab run change the seed to produce different random numbers

maxdE = 1;
diffLimit = 10^-5;
PrevE = zeros(1,N);
count = 0;

while (maxdE > diffLimit)
    
    v = randperm(N); %randon vecrot that have random values from 1 to N with no repetition, using it for fiding the energies with random permutations
%%%%%%%%%%%%%random update E vector with no consodiration for full sweep on the whole E vector once befor continueing to do second calculation on a specific site energy%%%%%%%%%%%%%
    for i=(1:N)
        E(v(i)) = epsilon(v(i)) + ((sum(r_nn*((1./((1+exp(E(1:v(i)-1)/T)))-0.5)./R_ij(v(i),1:v(i)-1))) + sum(r_nn*((1./((1+exp(E(v(i)+1:N)/T)))-0.5)./R_ij(v(i),v(i)+1:N)))));
    end
    
    if(count)
    DE = E - PrevE;
    maxdE = max(DE);
    end
    
    PrevE = E;
     
    count = count + 1;
    
end


