%Jonathan Braaten
%GM Project: Cation Exchange of Nafion 211
%2/20/2017-Final Modifications
%cation_exchange.m

clear;
clc;

%custom inputs
EW=1020; %gram/mol Equivalent weight of NR-211 Nafion (average of range 990-1050)
BW=50; %gram/m^2 Basic weight of NR-211 Nafion at 23 C and 50% RH
percH2O=5; %percentage of mass due to water (23C and 50%RH)-from NR-211 Spec Sheet (Avg between 2 and 8)


%material constants
t=25.4e-6; %m thickness of NR-211 Nafion at 23 C and 50% RH
w=.028; %m width of square sample
A_sample=w^2; %m^2 sample area 28 mm square
M_Co=291.03;    %g/mol molar mass of cobalt nitrate hexahydrate
M_H=63.01; %g/mol molar mass of nitric acid
p=.005; %Cation reduction factor from solution to membrane (0.1%)
Co_assay=0.98; %assay of cobalt nitrate hexahydrate (purity)
H_assay=0.70; %assay of nitric acid (purity)
rho_HNO3=1.413; %g/mL density of nitric acid at 20C 

%fit constants
a=0.0748;
b=2.513;
c=-1.041;
d=0.1976;

%calculated values
DW=BW*(1-(percH2O/100)); %Recalculating basic weight to account for water
m_sample=A_sample*DW; %grams mass of sample
N=m_sample/EW; %moles of sulfonic acid/hydrogen sites

%range of membrane charge fractions for exchange target
zeta_m=[0.6,0.7,0.8,0.9,0.92];

%calculate moles of material to add to solution
for i=1:length(zeta_m)
    zeta_s2(i)=(atan((zeta_m(i)-d)/a)-c)/b; %calculating equilibrium (s2) solution charge fraction
    M_H_m(i)=zeta_m(i)*N; %moles of protons in membrane after exchange
    M_Co_m(i)=((M_H_m(i)/zeta_m(i))-M_H_m(i))/2; %moles of cobalt ions in membrane after exchange 
    delta_M_Co(i)=M_Co_m(i); %moles of cobalt ions transferred to membrane during exchange process
    delta_M_H(i)=2*delta_M_Co(i); %moles of protons transferred to solution during exchange process (displaced by Cobalt ions)
    M_Co_s2(i)=(1/p)*delta_M_Co(i); %moles of cobalt ions in solution much greater than moles of cobalt ions transferred in exchange
    M_H_s2(i)=((2*zeta_s2(i)*M_Co_s2(i))/(1-zeta_s2(i))); %moles of protons in solution prior to exchange process
    M_Co_s1(i)= M_Co_s2(i)+delta_M_Co(i); %moles of cobalt ions in solution prior to exchange process
    M_H_s1(i)=M_H_s2(i)-delta_M_H(i); %moles of protons in solution after exchange process
    Mass_Co_Add(i)=(M_Co_s1(i)*M_Co)/Co_assay; %grams calculating mass of cobalt nitrate hexahydrate to add to solution
    Vol_HNO3_Add(i)=((M_H_s1(i)*M_H))/(rho_HNO3*H_assay); %mL calculating volume of nitric acid to add to solution
    zeta_s1_check(i)=M_H_s1(i)/(2*M_Co_s1(i)+M_H_s1(i)); %checking charge fraction for solution prior to exchange
    zeta_s2_check(i)=M_H_s2(i)/(2*M_Co_s2(i)+M_H_s2(i)); %checking charge fraction for solution after to exchange
    abs_error(i)=zeta_s1_check(i)-zeta_s2_check(i); %checking error between two zeta curves
    zeta_membrane_check(i)=(N-2*delta_M_Co(i))/(2*delta_M_Co(i)+(N-2*delta_M_Co(i))); %checking membrane charge fraction
end
figure(1)
plot(zeta_s2,zeta_m,':')
hold on
plot(zeta_s1_check,zeta_m,'x')
plot(zeta_s2_check,zeta_m,'+')
axis([0 1 0 1])
X=[0.0001,0.17,0.89,0.97,0.995,1];
Y=[0.07,0.145,0.36,0.74,0.86,0.95];
plot(X,Y,'o')
hold off
figure(2)
plot(zeta_m,zeta_membrane_check)
figure(3)
plot(zeta_m,abs_error)

