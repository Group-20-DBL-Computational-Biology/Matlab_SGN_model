clear all;
close all

%% input by user
%define timespan by user
tspan = 0:0.1:480;

%initial conditions
DNA0 = 3.4; %can be changed

%store the parameter values
para = [18.2 16.1 1.1e-2 7.8e-4 4.5e-3 0.2 8.5 65.8 6e-5];
    
%store the initial conditions
initCond = [DNA0 0  0 0 1 1];

%options for ODE solver
opts = odeset('RElTol', 1e-9, 'NonNegative', 1:6); 
%Explanation: NonNegative, because all six parameters cannot have a
%negative concentration. 
%% Solving the ODE system
ODEFUN = @(t,y) SystemState(t,y,para); %define your ODE function

[t_values, sol_values] = ode15s(ODEFUN, tspan, initCond, opts); 
%output will be in the following order: DNA dmRNAdt, dPdt, dPmatdt, dTsRdt,
%dTlRdt (DNA will be constant)

%convert minutes to hours
thours = t_values/60;

%% Plotting
%plot time vs  protein concentration and matured protein concentration
figure(1)
plot(thours, sol_values(:,3), thours,sol_values(:,4));
title('Concentration of protein in time');
xlabel('time (hours)');
ylabel('concentration (nM)'); 
legend('protein', 'matured protein')

%plot time vs mrna
figure(2)
plot(thours, sol_values(:,2));
title('Concentration of mRNA in time');
xlabel('time (hours)');
ylabel('concentration (nM)'); 

%plot time vs resources
figure(3)
plot(thours, sol_values(:,5), thours, sol_values(:,6));
title('Concentration of the resources in time');
xlabel('time (hours)');
ylabel('concentration (nM)'); 

%plot time vs GFP (NOT matured)
figure(4)
plot(thours, sol_values(:,3));
title('Concentration of protein in time');
xlabel('time (hours)');
ylabel('concentration (nM)'); 