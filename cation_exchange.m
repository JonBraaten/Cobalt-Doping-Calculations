%Jonathan Braaten
%GM Project: Cation Exchange of Nafion 211
%11/23/2016
%cation_exchange.m

clear;
clc;

%custom inputs
EW=1050; %gram/mol Equivalent weight of NR-211 Nafion
BW=50; %gram/m^2 Basic weight of NR-211 Nafion at 23 C and 50% RH
percH2O=0; % Percent percentage of mass due to water
V_H2O=.1; %L

%constants
t=25.4e-6; %m thickness of NR-211 Nafion at 23 C and 50% RH
w=.028; %m width of square sample
A_sample=w^2; %m^2 sample area 28 mm square
M_Co=291.03;    %g/mol molar mass of cobalt nitrate hexahydrate
M_H=63.01; %g/mol molar mass of nitric acid
AN=6.022e23; %part/mole Avogadro's Number
p=0.001; %acceptable reduction in concentration at equilibrium exchange (0.1%)
Co_assay=0.98; %assay of cobalt nitrate hexahydrate (purity)
H_assay=0.70; %assay of nitric acid (purity)
rho_HNO3=1.413; %g/mL density of nitric acid at 20C 

%calculated values
DW=BW/(1+(percH2O/100)); %Recalculating basic weight to account for water
m_sample=A_sample*DW; %grams mass of sample
N=m_sample/EW; %moles of sulfonic acid/hydrogen sites

%minimum and maximum zeta
Nco1=(0.08/2)*N;
Nh1=0.92*N;
Nco2=(0.3/2)*N;
Nh2=0.7*N;
zeta_m_max=Nh1/(2*Nco1+Nh1);
zeta_m_min=Nh2/(2*Nco2+Nh2);

zeta_m=zeta_m_min:.01:zeta_m_max; %zeta values for membrane

%Greszler curve fit parameters
a=0.0748;
b=2.513;
c=-1.041;
d=0.1976;

for i=1:length(zeta_m)
    zeta_s2(i)=(atan((zeta_m(i)-d)/a)-c)/b; %Greszler Curve fit to get solution's equilibrium charge fraction from the membrane's charge fraction 
    zeta_s1(i)=zeta_s2(i)/(1-p); %charge fraction prior to exchange
    n_H_m(i)=zeta_m(i)*N; %moles H+ in Nafion at equilibrium
    n_Co_m(i)=(N-n_H_m(i))/2; %moles Co2+ in Nafion at equilibrium
    n_Co_s(i)=(1/p)*n_Co_m(i); %moles Co2+ in solution (1/p times more than exchanged cobalt)
    n_H_s(i)=(2*zeta_s1(i)*n_Co_s(i))/(1-zeta_s1(i)); %moles H+ in solution
    m_Co_add(i)=(n_Co_s(i)/Co_assay)*M_Co; %grams Mass of cobalt nitrate hexahydrate to add to solution
    m_H_add(i)=(n_H_s(i)/H_assay)*M_H; %grams Mass of nitric acid to add to solution
    V_H_add(i)=m_H_add(i)/rho_HNO3; %mL volume of nitric acid to add to solution
    C_Co(i)=((m_Co_add(i)/M_Co)/V_H2O); %Mol/L Concentration of solution due to Cobalt cations
    C_H(i)=(m_H_add(i)/M_H)/V_H2O; %Mol/L concentration of solution due to protons
    zeta_check(i)=C_H(i)/(2*C_Co(i)+C_H(i)); %checking zeta of solution to reference with zeta_s1
end

plot(zeta_s1,zeta_m,'+:')

axis([0 1 0 1])
xlabel('\xi_{solution}')
ylabel('\xi_{membrane}')
hold on
X=[0.0001,0.17,0.89,0.97,0.995,1];
Y=[0.07,0.145,0.36,0.74,0.86,0.95];
plot(X,Y,'o')
legend('Experimental Prediction','Greszler Data Points')
figure(2)
plot(zeta_m,zeta_s1)
hold on
plot(zeta_m,zeta_check)

%this is a check for my branch
