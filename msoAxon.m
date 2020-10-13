
function[t,y]= msoAxon(stimType,start,stop,I,node,type,tEnd,v0, Input_node, Syn)
 
    switch(stimType)
        case{'Synaptic', 'SynapticPair'}
            [Syn.tSyn, Syn.gSyn] = Synaptic(Syn);
    end

    %%%Initial conditions
    m0 = .12; %Calculating with v0=-68 for all rest to get steady state vars
    h0 = .67; 
    p0 = 0; 
    w0 = .28; 
    z0 = .67; 
    a0 =.22;
    y0 = reshape(repmat([v0 m0 h0 p0 w0 z0 a0],45,1),45*7,1); 

    t = [];
    t0 = 0;
    y = []; 
    
    %%%Initial conditions
    load('Area');
    Syn.SA = SA; Syn.SAcm = SAcm; Syn.L = L; Syn.Lcm = Lcm; Syn.Rcyl = Rcyl; Syn.XA = XA; Syn.XAcm = XAcm;
    v0 =v0; 
    m0 = .12;
    h0 = .67; 
    p0 = 0; 
    w0 = .28; 
    z0 = .67; 
    a0 =.22;
    y0 = reshape(repmat([v0 m0 h0 p0 w0 z0 a0],45,1),45*7,1); 
    
    options = odeset('abstol',1e-8,'reltol',1e-8);
    %[t,y] =ode15s(@msoAxonODE, [0,tEnd], y0,[],stimType,start,stop,0,node,type,v0, Input_node, tEnd, Syn);
    %y0 = y(end,:);
    [t,y] =ode15s(@msoAxonODE, [0,tEnd], y0,[options],stimType,start,stop,I,node,type,v0, Input_node, tEnd, Syn);
end

%%%% ODEs FOR MODEL %%%%
function dx = msoAxonODE(t,x,stimType,start,stop,I,node,type,v0, Input_node, tEnd, Syn)

    % Inputs 
    % t: time [ms]
    % x: voltage and gating variables
    %
    % Outputs
    % dx: equations defining ODEs for voltage and gativng variables

    % Notes
    % equations as defined in Lehnert et al 2014 unless otherwise noted
    % all quantities expressed as dennsities, per unit area in sq cm [cm2]
    % 1 soma, 2 AIS, 21 internodes, 21 nodes = 45 compartments
    % dynamic variables: v, m+h [Na], p [KHT] , w+z[KLT] a [h]
    % possible to reduce: internodes are passive. KHT only in soma no h in nodes
    y = reshape(x,45,7);
    V = y(:,1); %0-45
    m = y(:,2); %46-90
    h = y(:,3); %91-135
    p = y(:,4); %136-180
    w = y(:,5); %181-225
    z = y(:,6);%226-270
    a = y(:,7); %271-315

    % specific capacitance. [multipy by 1e-8 to convert from uF/cm2 to uF/um2]  F = s/Ohm = s S. microF = ms mS, multiply by 1e6 to get ms nS
    C = [0.8 ; 0.8 ; 0.8 ; repmat([ 0.01 ; 0.8], 21, 1) ]*1e6*1e-8;  % from Lehnert, assuming that all unmyelinated areas have same capacitance
    SA = Syn.SA; SAcm = Syn.SAcm; L = Syn.L; Lcm = Syn.Lcm; Rcyl = Syn.Rcyl; XA = Syn.XA; XAcm = Syn.XAcm;
       
    % sodium [lehnert et al use Scott 2010]
    gNa = [0.2 ; 4 ; 4 ; repmat([0 ; 4], 21,1)]; % Na conductance, from Lehnert Table 2 [nS/um^2]
    VNa = 69;  % Na reversal potential [mV]
    INa  = gNa .* m.^4 .* (0.993*h +0.007) .* (V - VNa);  % Na current [pA/um^2]

    % high threshold potassium [lehnert et al use Rothman Manis 03 without slow activation]
    gKHT = zeros(size(V)); gKHT(1) = 0.1;  % KHT conductance only at soma, from Lehnert Table 2 [nS/um^2]
    VK = -90;  % K reversal potential [mV]
    IKHT =gKHT .* p .* (V - VK); % KHT current [pA/um^2]
    
    % low threshold potassium [lehnert et al use Mathews et al 2010]
    gKLT = [1.55 ; 1.55 ; 1.55 ; repmat([0 ; 1.55],21,1)]; % KLT conductance, from Lehnert Table 2 [nS/um^2]
    IKLT= gKLT .* w.^4 .* z .* (V - VK); % KLT current [pA/cm^2]

    % hyperpolarization activated cation current [lehnert et al use Baumann et al 2013]
    gh = [0.02 ; 0.02 ; 0.02 ; zeros(42,1)];  % h conductance, from Lehnert Table 2 [nS/um^2]
    Vh = -35; % h reversal potential [mV]
    Ih = gh .*  a .* (V - Vh); % h current [pA/cm^2] *** I don't see  power mentioned in Lehnert or Baumann so I make this linear in activation as in RM03 and others ***
    
    % leak current
    glk = [0.0005 ; 0.0005 ; 0.0005 ; repmat([0.0002 ; 0.05],21,1)];  % lk conductance, from Lehnert Table 2 [nS/um^2]
    Vlk = v0; % lk reversal potential [mV]  *** I don't see in Lehnert, so I set to stated resting potential of -68 mV
    Ilk = glk .* (V - Vlk); % lk current [pA/um^2] 

    % axial current 
    R = 100; % specific axial resistivity Ohm*cm. 
    Iaxial = zeros(size(V)); 
    Iaxial(1) = (2/R) * (V(2)-V(1)) / (Lcm(1)/XAcm(1) + Lcm(2)/XAcm(2)); % mV / (Ohm cm cm /cm^2) = mA
    Iaxial(2:44) = (2/R) * (V(1:43)-V(2:44)) ./ (Lcm(1:43)./XAcm(1:43) + Lcm(2:44)./XAcm(2:44)) + ...
                   (2/R) * (V(3:45)-V(2:44)) ./ (Lcm(3:45)./XAcm(3:45) + Lcm(2:44)./XAcm(2:44));
    Iaxial(45) = (2/R) * (V(44)-V(45)) / (Lcm(44)/XAcm(44) + Lcm(45)/XAcm(45)); 
    Iaxial = -Iaxial * 1e9 ./ SA; % convert mA to pA/um2
    % external current   
    switch(stimType)
        
        case('ramp')
            startRamp = start;    % time at start of ramp
            stopRamp = stop;     % slope of ramp [pA/ms]
            IrampMax  = I; % maximum current of ramp
            ramp= 1/(stopRamp-startRamp);
            tEnd = stopRamp + 5;
            tTop = startRamp+(1 / ramp);
            I0 = ((t>=startRamp)*(t<=tTop).*I.*(t-startRamp)*ramp)/1000;
            Iapp=-I0*1e3/SA(1);
            Iext = zeros(size(V)); 
            if t>5 && t<=tEnd; Iext(Input_node) =Iapp ; end
            
        case('ramp2')
            startRamp = start;    % time at start of ramp
            stopRamp = stop;     % slope of ramp [pA/ms]
            IrampMax  = I; % maximum current of ramp
            ramp2= 1/(stopRamp-startRamp);
            tEnd = tEnd;
            tTop = startRamp+(1 / ramp2);
            Ia = ((t>=startRamp)*I*(t-startRamp)*ramp2)/1000;
            Ib = I/1000;
            I0 = min(Ia,Ib);
            Iapp=-I0*1e3/SA(1);
            Iext = zeros(size(V)); 
            Iext(Input_node) =Iapp;
            
            
        case('step')
            Iext = zeros(size(V));
            if t>start && t<=stop; Iext(1) = -I/SA(1); end
            
         case('sine')
            Iext = zeros(size(V));
            I0 = @(T) (I *sin(2*pi*Syn.f*(T-start)/1000)* ...
           (sin(2*pi*Syn.f*(T-start)/1000)>0));
            Iapp = -I0(t)/SA(1); % convert to current density pA/um2
            if t>start && t<= stop; Iext(Input_node) = Iapp; end
            
         case('Synaptic')
            Iext = zeros(size(V)); % initialize
            g = interp1q(Syn.tSyn',Syn.gSyn',t);
            Iapp =  g * (V(1) - Syn.VsynE)/SA(1); 
            if t>start && t<= stop; Iext(Input_node) = Iapp; end
            
         case('SynapticPair')
            Iext = zeros(size(V)); % initialize
            g1 = interp1q(Syn.tSyn',Syn.gSyn',t); %Calculates Conductance
            g2 = interp1q(Syn.tSyn',Syn.gSyn',t+Syn.Diff);
            
            Iapp1 =  g1 * (V(1) - Syn.VsynE)/SA(1);
            Iapp2 =  g2 * (V(1) - Syn.VsynE)/SA(1);
            
            if t>start && t<= stop; Iext(Input_node) = Iapp1+Iapp2; end
            
        case('EPSG')
            Iext = zeros(size(V)); % initialize
            tEPSG = t-start;
            epsg = (1/0.21317) .* (tEPSG>0) .* (exp(-tEPSG/.18) - exp(-tEPSG/.1) ); % unitary epsg waveform
            Iapp = (I.*(0-V(Input_node)).*(tEPSG>=0).*epsg)/-SA(Input_node) ;
            if t>start && t<= stop; Iext(Input_node) = Iapp; end
        
        case('EPSGpair')
            Iext = zeros(size(V)); % initialize
            tEPSG = t-start;
            startEPSG1 = start;   % start of first EPSG event
            startEPSG2 = stop; % start of second EPSG event
            td = startEPSG2-startEPSG1;
            epsg = @(t,G) (G/0.21317) .* ones(size(t)) .* (t>0) .* (exp(-t/.18) - exp(-t/.1)); % unitary epsg waveform
            Iapp = (I.*(0-V(Input_node)).* ((tEPSG>=0).*epsg(tEPSG,1) + ...
                (tEPSG>=(td)).*epsg(tEPSG-td,1)))/-SA(Input_node);          
           Iext(Input_node)=Iapp;
    end
    
    % voltage dynamics %
    % all currents are density, unites pA/um2
    dV= - (INa + IKHT + IKLT + Ih + Ilk + Iext + Iaxial ) ./ C ; 
    
    % Na gating variable dynamics - Scott et al 2010 - 35 C
    minf = 1./(1 + exp((V+46)/(-11))); 
    taum = (  0.141 + (-0.0826 ./  (1 + exp((-20.5-V)/10.8)) )  ) / 3; 
    dm = (minf - m) ./ taum;
    
    hinf = 1 ./  (1 + exp((V+62.5)/7.77));
    tauh = (4 + ( -3.74 ./ (1 + exp((-40.6-V)/5.05))) ) / 3;
    dh = (hinf - h) ./ tauh;
    
    % KHT gating variables - Rothman Manis 2003 - inactivation only - 22 degrees
    % adjust to 35 C using Q10 of 3
    pTemp = 3^((22-35)/10);
    pinf = 1./ ( 1 + exp(-(V+23)/6));
    taup = pTemp * ( 100./(4*exp( (V+60)/32) + 5*exp( -(V+60)/22)) + 5);
    dp = (pinf - p) ./ taup;
    
    % KLT gating variables - Mathews et al 2010 - 35 degrees
    %TODO: copy+paste into EasyRun... don't delete, comment out
    winf = 1 ./ (1 + exp((V+57.34)/(-11.7)));
    tauw = 21.5 ./ ( 6*exp((V+60)/7) + 24*exp(-(V+60)/50.6)  ) + 0.35;
    dw = (winf - w) ./ tauw;
    
    zinf = (1-0.27)./(1 + exp( (V+67) / 6.16) ) + 0.27;
    tauz = 170 ./ ( 5*exp((V+60)/10) + exp(-(V+70)/8) ) + 10.7;
    dz = (zinf - z) ./ tauz;
    
    % h gating variables - Baumann et al 2013 - 32 degrees
    % adjust to 35 C using Q10 of 3
    aTemp = 3^((32-35)/10);
    ainf = 1./ (1+exp(0.1*(V+80.4)));
    taua = aTemp * ( 79 + 417*exp(-(V+61.5).^2/800));
    da = (ainf - a) ./ taua;
    
    % output
    type=type;
    switch(type)
        case('active-KLT')
            dx = [dV; 0*dm ; 0*dh ; 0*dp ; dw ; dz ; 0*da];
        case('active-H')
            dx = [dV; 0*dm ; 0*dh ; 0*dp ; 0*dw ; 0*dz ; da];
        case('active-KLT+H')
            dx = [dV; 0*dm ; 0*dh ; 0*dp ; dw ; dz ; da];
        case('passive')
            dx = [dV; dm*0 ; dh*0 ; dp*0 ; dw*0 ; dz*0 ; da*0]; %MAKES THIS LINEAR
        case('active-sodium')
            dx = [dV; dm ; dh ; 0*dp ; 0*dw ; 0*dz ; 0*da];
        case('active-full')
            dx = [dV; dm ; dh ; dp ; dw ; dz ; da];
        case('active-KHT')
            dx = [dV; dm*0 ; dh*0 ; dp ; dw*0 ; dz*0 ; da*0];
    end
end