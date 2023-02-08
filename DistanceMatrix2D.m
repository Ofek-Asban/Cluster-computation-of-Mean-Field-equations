function [R_ij, r_nn] = DistanceMatrix(N) %function that creats 2D sample with uniform distribured N sites and calculates the observed average nearest neighbors distance (r_nn)

%%%%%%%%%%%%%%%%%%% There is no need in this file for rng('shuffle') because SlaveHC2 have this and it calls DistanceMatrix.m so in different realizations there are different seed and thus different Matrices 

% creating random matrix for coordinats with nonrepeatable values
L=1*((N^(0.5))-1); %multiplying the size of the matrix by (N^0.5)-1 for making r_nn constant and thus making the density of sites rho = 1/(r_nn)^2 constant


c = L*(rand(N,2)-0.5); %c stand for coordinations, and multipling the system size by L=1

R_ij = zeros(N);

for i = 1:N
    for j = i+1:N
        
        dx = c(i,1)-c(j,1);
        dy = c(i,2)-c(j,2);
        
        %%%%%Implementing periodic boundary condistions, (PBC)%%%%%
        if abs(dx) > L/2
            dx = dx - L*sign(dx);
        end
        
        if abs(dy) > L/2
            dy = dy - L*sign(dy);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        distance = sqrt(dx^2+dy^2); %calculating the shortest euclidian distance between sites i and j 
        
        R_ij(i,j) = distance; % puting the distance in the distance matrix R_ij
        R_ij(j,i) = distance; % distance matrix always must be symmetric
        
    end
    
end

r_nn = 1; %calculating the analytical distance between nearest neighbors. which is constant (in this case r_nn=1)
%the observed r_nn in periodic boundary conditions is smaller in the same
%order of magnitude that the distance matrix is smaller in the PBC then the
%one with no PBC. Still it was observed numerically that r_nn in constant L
%and NO PBC is L/(N^0.5-1)

