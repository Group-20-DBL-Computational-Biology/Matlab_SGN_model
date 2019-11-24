%% Expression of a single protein
%This script models the expression of a single protein (GFP) from transcription
%to translation in-vitro. This model is based on the model made St?gbauer,
%Windhager, Zimmer and R?dler. Reactions can be found in the paper:
%Experiment and mathematical modeling of gene expression dynamics in a
%cell-free system.

%% Settings
timespan = 4;
number = 1000; %number of steps
dt = 1; %size of time step (minute)
if timespan ~= 0
    number = timespan/(dt/60);
end 
%Set the storage for the concentrations
mRNA = zeros(1,number);
P = zeros(1,number); %expressed protein
TsR = zeros(1,number); %transcription resource
TlR = zeros(1,number); %translation resource
Pmat = zeros(1,number); %maturated GFP (not necessary)

% Start concentrations (nanomolar)
mRNA(1) = 0; 
DNA = 3.4; %can be changed
P(1) = 0;
TsR(1) = 1; %can be changed
TlR(1) = 1; %can be changed
Pmat(1) = 0;

% Constants
kts = 18.2; %nanomolar per minute
ktl = 16.1; %nanomolar per minute
kcs = 1.1e-2; % per minute
deltamRNA = 7.8e-4; %per minute 
deltaTlr = 4.5e-3; %per minute
kmat = 0.2;%per minute
Ks = 8.5;% nanomolar
Kl = 65.8; %nanomolar
Ktlr = 6e-5; %

%array with the timesteps for plotting
time = 0:dt/60:(number)*dt/60; 
%% The Model
for i= 1:1:number
    %differential equations (see the paper)
    dmRNAdt = (kts * TsR(i) * DNA)/(Ks + DNA) - (deltamRNA * mRNA(i));
    dPdt = (ktl * TlR(i) * mRNA(i))/(Kl + mRNA(i)) - (kmat * P(i));
    dPmatdt = kmat * P(i);
    dTsRdt = -(kcs * TsR(i) * DNA)/(Ks + DNA);
    dTlRdt = -(deltaTlr * TlR(i)) /(Ktlr + TlR(i)); 
    
    %change of concentration for the timestep
    change_mRNA= dmRNAdt * dt; 
    change_P = dPdt *dt;
    change_Pmat = dPmatdt*dt;
    change_TsR = dTsRdt *dt;
    change_TlR =dTlRdt *dt; 
    
    %store new (total) concentration in arrays
    mRNA(i+1) = mRNA(i) + change_mRNA;
    P(i+1) = P(i) + change_P;
    Pmat(i+1) = Pmat(i) + change_Pmat;
    TsR(i+1) = TsR(i)+ change_TsR;
    TlR(i+1) = TlR(i) + change_TlR;
end 

%% Plotting the results
figure (1)
plot(time, P); 
title('Concentration of protein in time');
xlabel('time (hours)');
ylabel('protein (nM)'); 

figure(2)
plot(time, mRNA); 
title('Concentration of mRNA in time');
xlabel('time (hours)');
ylabel('mRNA (nM)'); 


