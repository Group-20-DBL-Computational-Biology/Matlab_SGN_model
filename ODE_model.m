close all
clear all;

%% input by user
%define timespan by user
tspan = 0:0.1:480;

%initial conditions
DNA0 = 3.4; %can be changed
L0 = 100e3; 
%store the parameter values
para = [18.2 16.1 1.1e-2 7.8e-4 4.5e-3 8.5 65.8 6e-5 8.33e-4 8.33e-4 60e13 0.8*60 60*7e18 60*6.2e-4 30e6*60 500*60 0.1*60];
    
%store the initial conditions
initCond = [DNA0 0 0 1 10 0 0 DNA0 0 L0 0 0 0];

%options for ODE solver
opts = odeset('RElTol', 1e-9, 'NonNegative', 1:13); 
%Explanation: NonNegative, because all six parameters cannot have a
%negative concentration. 
%% Solving the ODE system
ODEFUN = @(t,y) SystemState(t,y,para); %define your ODE function

[t_values, sol_values] = ode15s(ODEFUN, tspan, initCond, opts); 
%output will be in the following order: [dDNA; dmRNA_lacI; dmRNA_lacZ; dTsR; dTlR; dRA; dR; dO; dRO; dL; dG; dGL; dA];


%convert minutes to hours
thours = t_values/60;
%% Plotting
%plot time vs  protein concentration
figure(1)
plot(thours, sol_values(:,7), thours,sol_values(:,11));
title('Concentration of protein in time');
xlabel('time (hours)');
ylabel('concentration (nM)'); 
legend('repressor', 'galactosidase')

%plot time vs mrna
figure(2)
plot(thours, sol_values(:,2),thours, sol_values(:,3) );
title('Concentration of mRNA in time');
xlabel('time (hours)');
ylabel('concentration (nM)');
legend('mRNA lacI', 'mRNA lacZ')

%plot time vs resources
figure(3)
plot(thours, sol_values(:,4), thours, sol_values(:,5));
title('Concentration of the resources in time');
xlabel('time (hours)');
ylabel('concentration (nM)'); 
legend('TsR', 'TlR')

figure(4)
plot(thours, sol_values(:,10), thours, sol_values(:,13));
title('Concentration of the resources in time');
xlabel('time (hours)');
ylabel('concentration (nM)'); 
legend('Lactose', 'Allolactose')
%%
figure(5)
plot(thours, sol_values(:,7), thours, sol_values(:,6));
title('Concentration of the resources in time');
xlabel('time (hours)');
ylabel('concentration (nM)'); 
legend('Repressor', 'complex RA')


