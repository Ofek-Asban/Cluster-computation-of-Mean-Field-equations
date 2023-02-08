function [] = HcTLS(RealNum, N) % Gets a number of the specific realization num (RealNum) adn the system size N, heating and cooling the sample between T min to T max in 20 steps, and then saves the solution E vector into a file with the specific realization number

rng('shuffle') %This command in each matlab run change the seed to produce different random numbers (Only relevant on distributed computing because each jobs that goes to new porcessor opens matlab for the first time which gives same seed that generetes same random numbers)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Temperature%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T_low = 1;
%initial temperature which will be the final temperature at the end of the process (0.05 in Ariel's paper for calculation of the DOS and 0.1 for calculation of the distribution of eigen values)
T_high = 10; %hottest temperature (in this case 1 bcs there the temp will be comperable with the interactions though u_ij/(T*(r_nn)^3)=20)
%RealNum - number of specific realization

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Disorder (W=std of delta)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delta = normrnd(0,10,1,N); %1st num=mean 2nd num=std (disorder W), 3rd and 4th num = matrix dimensions

%calibrating "random on-site energies" (delta vector is the TLSs Asymmetry energy without the interactions-like epsilon_i in the electrin glass case only in this case the distribution is Gaussian and not uniform and it souldnt change much). 
Delta = 20*(rand(1,N)-0.5); % initial mean-field "average" TLS assymetry energy vertor for input

cd /fastspace/users/asban/parallelTLS %Changing the directory in fastspace for loading the Dmatrix.mat
[R_ij, r_nn] = DistanceMatrix3D(N); %creating in each realization new distance matrix

%%%%%%%%%%%Interactions (J=std of C_ij (where C_ij is u_ij in litrature))%%%%%%%%%%%
C_ij = normrnd(0,1,N,N); 
%creating a frefactor (in the shape of NxN matrix because of vectorization method in matlab) that the interaction term will have. With mean 0 and std 1. this term effect's the DOS (moshes paper 1307.0868v3) quntitevily, in case of uniform distribution for example we get the same gap but slightly different values.
C_ij = C_ij - tril(C_ij,-1) + triu(C_ij,1)'; %creating the matrix symmetric as it should be. for example bonds C_ij(3,1) = C_ij(1,3), as with distance matrix diagonal elements dont have any meaning (in Dmatrix they are just 0 here we have some unused value). 

p=10^(-3); %p is from p/delta_0 distribution of delta_0 (p helps control the cutoff range of delta_0 variable). According to esquinazi in all amorphous materials p*<|C_ij|> ~ 10^(-3) while J_ij=C_ij/(r_ij)^3 and <|J_nn|> = <|C_nn|>/(r_nn)^3 = <|C_ij|>/(r_nn)^3 = 1 we set r_nn = 1 so we get <|C_ij|> = 1 therefore p=10^(-3) to 10^-4

%making a random vector size N that represents delta_0 (tuneling energy) with 1/delta_0 distribution
n=10^(p*6);%the base of the log which together with 'a' will give the range of the distribution of 1/delta_0.
a=10;

%The range of the random vector of uniform distribution that with it I create the vector with 1/delta_0 
y1=(log(1/(a*(n^(1/p)))))/(log(n));
y2=(log(1/a))/(log(n));

%building the tunneling energy vector with 1/delta_0 distribution out of a vector that is uniform distributed in the range [1/a*n^(1/p), 1/a]
v = (abs(y2-y1))*rand(1,N)+min(y1,y2); 
delta_0 = n.^(v);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%Heating and then Cooling procedure%%%%%%%%%%%%%%%%%%%%%%%%%% 
T = T_low; % initializing the temperature
Num_Of_T_Steps = 10; %number of heating and thien cooling steps (ref 18 used 20 cooling steps)
size_Of_T_Step = 2*(T_high-T_low)/(Num_Of_T_Steps); %size of each step (multiplyed by 2 because we need to heat and then cool down the system)

for k = 1:(Num_Of_T_Steps+1)
    
    [Delta] = EnergyIterations(N, delta, Delta, delta_0, R_ij, T, r_nn, C_ij); %calling the function to get the solution of E vector for a given T and saving it in a file
     
    if(k <= Num_Of_T_Steps/2) %Heating the system
        T = T + size_Of_T_Step;
    end
    
    if(k > Num_Of_T_Steps/2) %Cooling back the system
        T = T - size_Of_T_Step;
    end	
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%Getting the final energy E of the TLS -> LOOK Tunneling Systems in Amorphous and Crystalline Solids Pablo Esquinazi PAGES 234-243 & 279-280%%%%%%%%%%%%%%%
    E = sign(Delta).*sqrt(Delta.^2 + delta_0.^2); %finding the energy difference of the diagonal TLS after finding the Delta of the TLS and taking the unchaned tunneling energy Delta_0 (Delta_0 is static and plays no role in the interations because we know the in a good approx Delta_0 doesnt change under the interaction with the strain field)

cd /fastspace/users/asban/output_mat %changing again the directory to output_mat to save there the SleaveHeatingAndCooling_RealNum file
fname = sprintf('Eanddelta0_%d.mat', RealNum); %adding the realization number to the file name
save(fname, 'E', 'delta_0'); %saving the E vector to a new file of current Realization (these E vector in the file will represent a realization of the system after heating and cooling the system).
%fname = sprintf('Edelta0DmatWJ_%d.mat', RealNum); %adding the realization number to the file name
%save(fname, 'E', 'delta_0', 'Delta', 'r_nn', 'R_ij', 'C_ij'); Saving also R_ij and C_ij in order to calculate the full rate equation with the direct interaction term
