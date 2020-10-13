function [t,x]=TwoCpt(stimType,start,stop,I,node,type,tEnd,v0,inputNode, Syn)
load('Area.mat');
P = getParam(v0, node, inputNode, type);

%%% Run two-compartment model 
switch(stimType)
    case('step')
        startStep = start;    % time at start of step
        stopStep  =  stop;  % time at end of step
        Istep     = I; % current level of step 
        P.step =Istep; P.startStep = startStep; P.stopStep = stopStep; odeFile = @TwoCptODE;
    case('ramp')
        startRamp = start;    % time at start of ramp
        stopRamp = stop;     % slope of ramp [pA/ms]
        IrampMax  = I; % maximum current of ramp
        P.startRamp =startRamp; P.ramp= 1/(stopRamp-startRamp); P.I = IrampMax; odeFile = @TwoCptODE;
    case('ramp2')
        startRamp = start;    % time at start of ramp
        stopRamp = stop;     % slope of ramp [pA/ms]
        IrampMax  = I; % maximum current of ramp
        P.startRamp =startRamp; P.ramp2= 1/(stopRamp-startRamp); P.I = IrampMax; odeFile = @TwoCptODE;
    case('sine')
        Ieq  = I; % Equation of current
        startSine = start;    % time at start of step
        stopSine =  stop;  % time at end of step
        P.sine =Ieq; P.startSine = startSine; P.stopSine = stopSine; P.f = Syn.f; odeFile = @TwoCptODE;
    case('Synaptic')
        [Syn.tSyn, Syn.gSyn] = Synaptic(Syn);% Calculates
        P.tSyn=Syn.tSyn; P.gSyn=Syn.gSyn; P.stopSyn=stop; P.startSyn=start; 
        P.VsynE=Syn.VsynE; odeFile = @TwoCptODE;
    case('SynapticPair')
        [Syn.tSyn, Syn.gSyn] = Synaptic(Syn);% Calculates
        P.tSyn=Syn.tSyn; P.gSyn=Syn.gSyn; P.stopSyn2=stop; P.startSyn=start; 
        P.VsynE=Syn.VsynE; P.Diff=Syn.Diff; P.tEnd = tEnd; odeFile = @TwoCptODE;
    case('EPSG')
        startEPSG = start;    % start of EPSG event
        gEPSG     =  I;   % amplitude of EPSG event
        P.EPSG=1;P.startEPSG =startEPSG; P.I= gEPSG; odeFile = @TwoCptODE;
    case('EPSGpair')
        startEPSG1 = start;   % start of first EPSG event
        startEPSG2 = stop; % start of second EPSG event
        gEPSG      = I;  % amplitude of EPSG event
        P.EPSGpair=1;P.startEPSG1 =startEPSG1; P.td = startEPSG2-startEPSG1; P.I= gEPSG; odeFile = @TwoCptODE;
end

options = odeset('abstol',1e-6,'reltol',1e-6,'maxstep',.1);
[t,x] =ode15s(odeFile, [0 tEnd], [P.Vrest P.Vrest P.winf(P.Vrest) ...
    P.hinf(P.Vrest) P.winf(P.Vrest) P.minf(P.Vrest) P.minf(P.Vrest) ...
    P.hinf(P.Vrest) P.pinf(P.Vrest) P.ainf(P.Vrest) P.ainf(P.Vrest)],options,P);
end

function P = getParam(v0, node, Input_node, type)
%%% Loading Constants
load('Coupling.mat');
load('Fractions.mat');
load('Area.mat');

%%% Finds Coupling Constants
P.couple12 = coupling1(node-1); % Forward coupling
P.couple21 = coupling2(node-1); % Backward coupling

%%% Calculates KLT Frac
switch (type)
    case{'active-KLT', 'active-full', 'active-KLT+H'}
        KLTfrac1 = KLT_frac(1); % KLT in compartment 1
        KLTfrac2 =KLT_frac(node); % KLT in compartment 2
    otherwise
        KLTfrac1 = 0; % KLT in compartment 1
        KLTfrac2 = 0; % KLT in compartment 2     
end

%%% Calculates H Frac
switch (type)
    case{'active-H', 'active-full', 'active-KLT+H'}
        hfrac1 = H_frac(1); % H in compartment 1
        hfrac2 = H_frac(node); % H in compartment 2
    otherwise
        hfrac1 = 0; % H in compartment 1
        hfrac2 = 0; % H in compartment 2     
end

%%% Calculates Na Frac
switch (type)
    case{'Active-sodium', 'active-full'}  
        NaFrac1 =  Na_frac(1); % Na in compartment 1
        if node == 3
            P.gNa2= 119;  % Na in compartment 2
        elseif node == 5
            P.gNa2= 25.5; % Na in compartment 2
        else
            P.gNa2= 140; %This Value is a placeholder for non compartment 2 Na
        end 
    otherwise
        NaFrac1 = 0; % Na in compartment 1
        P.gNa2=0;    % Na in compartment 2 
end

%%% Calculates KHT
switch (type)
    case{'active-KHT', 'active-full'}
        P.gKHT = 0.1 /SA(1) * 1000;
    otherwise
        P.gKHT = 0;
end

%%%% Fixed parameters %    HARD CODED PARAMETERS
P.areaRatio = areaRatio(node); % CPT1 to CPT2 area ratio
P.R1    = 10 * 1e-3;  %8.5  % Input resistance to CPT1 [10^9 Ohm]
tauEst  = .71;   % "estimated time constant" [ms]
Vrest   = v0;   %-58 % resting membrane potential [mV]
Elk     = v0;   %-58 % leak reversal potential [mV]
Vk = -90; % resting membrane potential for KHT

% Passive parameters %
gC    = (1/P.R1) * P.couple21 / (1-P.couple12*P.couple21); % coupling conductance [nS]
gTot1 = gC * (1/P.couple21 - 1); % CPT1 total conductance [nS]
gTot2 = gC * (1/P.couple12 - 1); % CPT2 total conductance [nS]

% Passive parameters that require separation of time scales assumption %
tau1  = tauEst * (1-P.couple12*P.couple21);          % CPT1 time constant [ms]
tau2  = tau1   * P.areaRatio * (P.couple12/P.couple21);  % CPT2 time constant [ms]
cap1  = tau1   * (gTot1 + gC); % CPT1 capacitance [pF]
cap2  = tau2   * (gTot2 + gC); % CPT2 capacitance [pF]

%%%%% PARAMETER STRUCTURE %%%%
P.areaRatio= areaRatio; P.tauEst= tauEst;   P.tau1  = tau1; P.tau2  = tau2;
P.cap1     = cap1; P.cap2  = cap2;     P.gTot1 = gTot1;    P.gTot2 = gTot2;
P.gC       = gC;   P.Vrest = Vrest;    P.Elk   = Elk;  P.VK    = Vk;

%%% KLT CONDUCTANCE %%%
P.EK   = -90;
P.winf = @(V)  (1. + exp(-(V+57.34)/11.7) )^-1; 
P.zinf = @(V) (1-0.27)./(1 + exp( (V+67) / 6.16) ) + 0.27;
P.tauw = @(V) 21.5 ./ ( 6*exp((V+60)/7) + 24*exp(-(V+60)/50.6)  ) + 0.35;
P.tauz = @(V) (170 ./ ( 5*exp((V+60)/10) + exp(-(V+70)/8) ) + 10.7);

P.gKLT1 = (KLTfrac1*P.gTot1) / (P.winf(P.Vrest)^4*P.zinf(Vrest));
P.gKLT2 = (KLTfrac2*P.gTot2) / (P.winf(P.Vrest)^4*P.zinf(Vrest));

%%% H Conductance %%%
P.Vh = -35; % h reversal potential [mV]

P.aTemp = 3^((32-35)/10);
P.ainf = @(V) 1./ (1+exp(0.1*(V+80.4)));
P.taua = @(V) P.aTemp * ( 79 + 417*exp(-(V+61.5).^2/800));

P.gh1 = (hfrac1 * P.gTot1) / P.ainf(P.Vrest); 
P.gh2 = (hfrac2 * P.gTot2) / P.ainf(P.Vrest); 

%%% Calculate Leak Currents %%%
P.glk1 = (1-KLTfrac1 - hfrac1- NaFrac1)*P.gTot1;
P.glk2 = (1-KLTfrac2 - hfrac2)*P.gTot2;

%KHT Voltage
P.pTemp = 3^((22-35)/10);
P.pinf =@(V) 1./ ( 1 + exp(-(V+23)/6));
P.taup =@(V) P.pTemp * ( 100./(4*exp( (V+60)/32) + 5*exp( -(V+60)/22)) + 5);

% Na gating Rothman Manis with 35C temp adjustment
P.taum = @(V) ((0.141 + (-0.0826 ./  (1 + exp((-20.5-V)/10.8)) )  ) / 3); 
P.minf = @(V)(1.+exp((V+46.)/-11.)).^-1.;
P.hinf = @(V)(1.+exp((V+(62.5))/7.77)).^-1.; 
P.tauh = @(V) (4 + ( -3.74 ./ (1 + exp((-40.6-V)/5.05))) ) / 3;

P.gNa1 = (NaFrac1 * P.gTot1) / (P.minf(P.Vrest)^4 * (0.993*P.hinf(P.Vrest) +0.007)) ; 
P.ENa = 69;     % Na reversal potential [mV]
P.gSyn = 26.7; % 6mV EPSP in V1 with Franken EPSG waveform

%Changes Input of Current Dependent on Input Location
if Input_node == 1
    P.IappLoc = 1;
else
    P.IappLoc = 2;
end
end