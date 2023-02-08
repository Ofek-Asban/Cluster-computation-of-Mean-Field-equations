function [Delta] = EnergyIterations(N, delta, Delta, delta_0, R_ij, T, r_nn, C_ij) %this function number of sites N calling the function distanceMatrixandRnn and calculated metastable solution to eq(5) and returns the solution in a shape of a vectors E for average site enegies and n vector of site occupations

maxdDelta = 1;
diffLimit = 10^-4;
PrevDelta = zeros(1,N);
count = 0;

while (maxdDelta > diffLimit)
    
    v = randperm(N); %randon vector that have random values from 1 to N with no repetition, using it for fiding the energies with random permutations
    %%%%%%%%%%%%%random update Delta vector with no considiration for full sweep on the whole Delta vector once befor continueing to do second calculation on a specific site energy%%%%%%%%%%%%%
    for i=(1:N)
        
        %E = sign(Delta).*sqrt(Delta.^2 + delta_0.^2); %passing on the sign of delta to E because E is 1/2(E_1-E2) sometimes E_1 > E2 so 1/2(E_1-E2)=E and this is when Delta > 0 and sometimes E_1 < E2 so 1/2(E_1-E2)= -E and this is when Delta < 0 
        %Delta(v(i)) = delta(v(i)) + sum((r_nn^3)*((C_ij(v(i),1:v(i)-1)).*((1./((1+exp(E(1:v(i)-1)/T)))-0.5)./((R_ij(v(i),1:v(i)-1)).^3)))) + sum((r_nn.^3)*((C_ij(v(i),v(i)+1:N)).*((1./((1+exp(E(v(i)+1:N)/T)))-0.5)./((R_ij(v(i),v(i)+1:N)).^3))));
        Delta(v(i)) = delta(v(i)) + 2*sum((r_nn^3)*((C_ij(v(i),1:v(i)-1)).*((1./((1+exp((sign(Delta(1:v(i)-1)).*sqrt(Delta(1:v(i)-1).^2 + delta_0(1:v(i)-1).^2))/T)))-0.5)./((R_ij(v(i),1:v(i)-1)).^3)))) + 2*sum((r_nn.^3)*((C_ij(v(i),v(i)+1:N)).*((1./((1+exp((sign(Delta(v(i)+1:N)).*sqrt(Delta(v(i)+1:N).^2 + delta_0(v(i)+1:N).^2))/T)))-0.5)./((R_ij(v(i),v(i)+1:N)).^3))));
    end
        
    if(count)
        DDelta = Delta - PrevDelta;
        maxdDelta = max(DDelta);
    end
    
    PrevDelta = Delta;
    
    count = count + 1;
    
end
