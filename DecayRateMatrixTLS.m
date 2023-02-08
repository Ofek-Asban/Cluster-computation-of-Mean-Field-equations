%Rates od a specific TLS system realization INCLUDING THE INTERACTIONS 
%This func will take an enrgy vector of a specific realization and (1) calculate the transition matrix A_ij WITH THE INTERACTION TERM of the master equation
%and diagonalize it and create a file with the eigen values vector "lambda"
%Where E is the energy vetor, T temperature which in our case for N=1000 should be 0.1 (like in the paper for finding the eigenvalues distribution), RealNum- realization number so we could create the suitible file name

N=1000;
T = 0.1;
J_nn = 0.3; %average near neighbor interaction strangth in Kelvins

%cd /media/asban/CDROM/Edelta0Dmat %changing the directory so that the for loop could reach all the slaveheatingandcooling files and for saving the result to a new file.
%filename = ['DmatEDelt0DeltPrime_', num2str(RealNum)]; % adding to the the filename 'SlaveHeatingAndCooling_' a string of the corrent index i of the loop to get the to i'th realization file
%Data = load(filename, 'E', 'delta_0'); % loading struct that contains the file data (R_ij, r_nn, E)                      
cd /users/asban/desktop/
Data = load('WJ', 'E', 'delta_0', 'Delta', 'r_nn', 'R_ij', 'C_ij');

E = [Data().E]; % Extracting the vector E from the struct E()
delta_0 = [Data().delta_0]; %Tunnneling energy
Delta = [Data().Delta]; %Meanfield assymetry energy (usually marked as Delta')
u = [Data().C_ij]; %Interaction prefactor (which holds the interaction strangth i.e. sqrt(<C_ij^2>) = U_0)
r = [Data().R_ij]; %distance matrix

%ALWASY CHECK: that r_nn is set to 1 when creating a distance matrix
%otherwise it will be needed to normalize the interaction term accorsing to
%the new J_nn= U_0/(r_nn)^3

A = zeros(N,N);
% Building the transition matrix A_ij for i!=j (only off diagonal elements). i and j running on the upper triangle not including the diagonal, and inputing on the upper triangle and down triagle at the same iteration
for i = 1:N %runing on rows (TLS i)
    C_i = (10^8)*(J_nn^3)*((delta_0(i))^2)*(Delta(i))*((sinh(E(i)/T))^(-1)); %Some prefactor that turns the coefficiants of A matrix to equality to units of J_nn (i.e a_i*J_nn)
    for j = 1:N %running on columns, for all the other TLSs that interact with TLS i (only the off diagonal terms)
        if(j==i)
            continue;
        end
        A(i,j) = (C_i*(u(i,j))/(T*(r(i,j))^3)); 
    end       
end

%After finishing bouilding the off-diagonal elements of A_ij, we are building the diagonal elements, making sure that current is preserved i.e. the sum of each line in A_ij is zero.
for i= 1:N
    A(i,i) = -(10^8)*(J_nn^3)*((delta_0(i))^2)*abs(E(i))*coth(abs(E(i))/(2*T)); %Each TLS rates lambda_i
end

lambda_full = real(eig(A)); %calculaing and saving in vector lambda which is the real part of the eigenvalues of FULL TLS rate matrix A including the direct interactions tern of the rate equation (the imaginary terms ,if we look on the decoupled rate eaquations, will give oscilatory terms inside the envelope of the decaying exponent).
lambda_diag = -(10^8)*(J_nn^3)*((delta_0).^2).*abs(E).*coth(abs(E)./(2*T)); %only the diagonal of the rate matrix which doesnt contain the direct interaction terms in the rate matrix (contains interactions only though the MF energies vector E)
lambda_full = sort(lambda_full);
lambda_diag = lambda_diag';
lambda_diag = sort(lambda_diag);

Rel_err = abs((lambda_full-lambda_diag)./lambda_full); %Relative error between the rates only on the diagonal to the full rate including the oof diagonal element (i.e direct interaction terms in the rate equation)

x=1:1:1000;
y=0.1*ones(N);
%plot(x,lambda_diag, '*',x,lambda_full, '+');
plot(x, Rel_err, '*', x, y)
