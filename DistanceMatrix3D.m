function [R_ij, r_nn] = DistanceMatrix3D(N) %function that creats 3D sample with uniform distribured N sites with constant site density.

% creating random matrix for coordinats with nonrepeatable values
L=1*((N^(1/3))-1);%In 3D lattice multiplying the size of the matrix by (N^(1/3))-1 for making r_nn constant and thus making the density of sites rho = 1/(r_nn)^3 constant aswell

c = L*(rand(N,3)-0.5); %c stand for coordinations, and multipling the system size by L, -0.5 comes from the fact that the cube will stance simatrically about the origen

R_ij = zeros(N);

for i = 1:N-1
    for j = i+1:N
        dx = c(i,1)-c(j,1); %index is site number and number is coordinate number (c(i,1) - x coordinate for site i, c(j,1) - x coordinate for site j)
        dy = c(i,2)-c(j,2);
        dz = c(i,3)-c(j,3);
            
        %%%%%Implementing periodic boundary condistions, (PBC)%%%%%
        if abs(dx) > L/2
            dx = dx - L*sign(dx);
        end
        
        if abs(dy) > L/2
            dy = dy - L*sign(dy);
        end
        
        if abs(dz) > L/2
            dz = dz - L*sign(dz);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        distance = sqrt(dx^2+dy^2+dz^2); %calculating the shortest euclidian distance between sites i and j 
        
        R_ij(i,j) = distance; % puting the distance in the distance matrix R_ij
        R_ij(j,i) = distance; % distance matrix always must be symmetric
        
    end
    
end


r_nn = L/((N^(1/3))-1); %calculating the analytical distance between nearest neighbors. which is constant (in this case r_nn=1 which makes a constant density of sites in the lattice, in 3D rho = 1/(r_nn^3))
        
