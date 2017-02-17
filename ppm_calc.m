%Jonathan Braaten
%GM Project: ICP-MS Capability Check
%10/26/2016
%ppm_calc.m

clear;
clc;

t1=25.4e-6; %m thickness of NR-211 Nafion at 23 C and 50% RH
ts=1000e-9; %m thickness of NR-211 spun coat on polycarb (thinnest layer)
ds=10e-3; %m diameter of NR-211 spun coat on polycarb 
Vs=pi*(ds/2)^2*ts; %m^3 volume of NR-211 sample
EW=1050; %gram/mol Equivalent weight of NR-211 Nafion
BW=50; %gram/m^2 Basic weight of NR-211 Nafion at 23 C and 50% RH
dens=BW/t1; %gram/m^3 Basic density of NR-211 Nafion at 23 C and 50% RH
percH2O=8; % Percent percentage of mass due to water
DW=BW/(1+(percH2O/100)); %Recalculating basic weight to account for water
A_sample=0.028*0.028; %m^2 sample area 28 mm square
m_sample=dens*Vs; %grams mass of sample
M=m_sample/EW; %moles of sulfonic acid/hydrogen sites
M2=M/2; %moles of cation sites (2 cobalt ions per sulfonate group)
M_co=58.933; %g/mol molar mass of cobalt
Fill=[1,10,50,100]; % percentage of sulfonate sites occupied by cobalt cations

for i=1:length(Fill)
    mass_cobalt(i)=((M2*Fill(i)/100)*M_co)*1000; %mass Cobalt in milligrams
    ppm15(i)=mass_cobalt(i)/0.015; %concentration ppm considered in 15 mL of water
    ppb15(i)=ppm15(i)*1000; %concentration ppb considered in 15 mL of water
    ppt15(i)=ppb15(i)*1000; %concentration ppm considered in 15 mL of water
end


