function dx = TwoCptODE(t,x,P)

    % Inputs %
    % t = time [ms]
    % x = Voltage and gating variables
    % P = structure of parameter values
        
    %%%% votage variables %%%% 
	V1  = x(1); 
    V2  = x(2);  
    
    %%%%% gating variables %%%% 
    w1  = x(3);
    z1  = P.zinf(P.Vrest);
    m1   = x(6);%P.minf(V1);
    m2   = x(7);%P.minf(V2);
    h1   = x(4);
    w2  = x(5);
    z2  = P.zinf(P.Vrest);
    h2   = x(8);
    p   = x(9);
    a1  = x(10);
    a2  = x(11);

    %%%%% INPUT CURRENT TO CPT 1 %%%% 
    if isfield(P,'step')

        if (t>=P.startStep && t<P.stopStep); Iapp = P.step; else Iapp = 0; end; 
        
    elseif isfield(P,'ramp')
       if ~isfield(P,'startRamp'); P.startRamp = 5; end;
       tTop = P.startRamp+(1 / P.ramp);
       Iapp = (t>=P.startRamp)*(t<=tTop).*P.I.*(t-P.startRamp)*P.ramp;
       
    elseif isfield(P,'EPSG')
        tEPSG = P.startEPSG; 
        epsg = @(t,G) (G/0.21317) .* ones(size(t)) .* (t>0) .* (exp(-t/.18) - exp(-t/.1) ); % unitary epsg waveform
        Iapp = P.I.*(0-V1).*(t>=tEPSG).*epsg(t-tEPSG,1) ;
        
    elseif isfield(P,'EPSGpair')
        tEPSG = t-P.startEPSG1;   
        epsg = @(t,G) (G/0.21317) .* ones(size(t)) .* (t>0) .* (exp(-t/.18) - exp(-t/.1) ); % unitary epsg waveform
        
        Iapp = (P.I.*(0-V1).* ((tEPSG>=0).*epsg(tEPSG,1) + ...
                (tEPSG>=(P.td)).*epsg(tEPSG-P.td,1)));             
           
    elseif isfield(P,'sine')
       Iext = @(T) (P.sine*sin(2*pi*P.f*(T-P.startSine)/1000)* ...
           (sin(2*pi*P.f*(T-P.startSine)/1000)>0));
       Iapp = (t>=P.startSine) *  (t<=P.stopSine) * Iext(t);
       
    elseif isfield(P,'stopSyn')
       g = interp1q(P.tSyn',P.gSyn',t);
       Iext =  g * (V1 - P.VsynE );
       Iapp = -(t>=P.startSyn) *  (t<=P.stopSyn) * Iext;
       
    elseif isfield(P,'stopSyn2')
       g1 = interp1q(P.tSyn',P.gSyn',t);
       if t+P.Diff < P.tEnd
        g2 = interp1q(P.tSyn',P.gSyn',t+P.Diff);
       else
        g2 = 0;
       end

       Iext1 =  g1 * (V1 - P.VsynE);
       Iext2 =  g2 * (V1 - P.VsynE);
       Iapp = -(t>=P.startSyn) *  (t<=P.stopSyn2) * (Iext1 + Iext2);           
    elseif isfield(P,'ramp2')
       P.startRamp = 5;       
       Iapp_2 = (t>=P.startRamp)*P.I.*(t-P.startRamp)*P.ramp2;
       Iapp = min(P.I,Iapp_2);

    end
     %%%%% Cpt1 currents %%%%% 
    IKHT  = P.gKHT .* p .* (V1 - P.VK);
    INa0  = P.gNa1 .* P.minf(P.Vrest)^4 * (0.993*P.hinf(P.Vrest) +0.007) .* (P.Vrest-P.ENa); %change into 1 eq
    INa1  = P.gNa1 .* m1^4 * (0.993*h1 +0.007) .* (V1-P.ENa) -INa0;
    Ilk1  = P.glk1 * (V1 - P.Elk);
    IKLT0 = P.gKLT1 * P.winf(P.Vrest)^4*P.zinf(P.Vrest)*(P.Vrest-P.EK);
    IKLT1 = P.gKLT1 * w1^4*z1*(V1-P.EK) - IKLT0;
    Ih0   = P.gh1 .*  P.ainf(P.Vrest) * (P.Vrest - P.Vh);
    Ih1   = P.gh1 .*  a1 .* (V1 - P.Vh) -Ih0;

    %%%%% Cpt2 currents %%%%%
    Ilk2  = P.glk2  * (V2 - P.Elk);
    INa0  = P.gNa2 .* P.minf(P.Vrest)^4 * (0.993*P.hinf(P.Vrest) +0.007) .* (P.Vrest-P.ENa); %change into 1 eq
    INa2  = P.gNa2 .* m2^4 * (0.993*h2 +0.007) .* (V2 -P.ENa) -INa0;
    IKLT0 = P.gKLT2 * P.winf(P.Vrest)^4*P.zinf(P.Vrest)*(P.Vrest-P.EK);
    IKLT2 = P.gKLT2 * w2^4*z2*(V2-P.EK) - IKLT0;
    Ih0   = P.gh2 .*  P.ainf(P.Vrest) .* (P.Vrest - P.Vh);    
    Ih2   = P.gh2 .*  a2 .* (V2 - P.Vh) -Ih0;
    
    %%%%% Coupling current %%%%%
    IC = P.gC*(V1-V2);
    
    %%%%% Update V using Current Balance Equation %%%%%
    
    if P.IappLoc == 1 %% Changes input to Soma or Axon
        Iapp1 = Iapp;
        Iapp2 = 0;
    else
        Iapp1 = 0;
        Iapp2 = Iapp;
    end
        
    dV1 =  ( -Ilk1 - IKLT1 - IC + Iapp1 - INa1 -Ih1 -IKHT)/ P.cap1 ;
    dV2 =  ( -Ilk2 - IKLT2 + IC + Iapp2 - INa2 -Ih2  )/ P.cap2 ; 
     
    %%%%% Update Gating Vars %%%%%
    dw1  = (P.winf(V1) -w1) / P.tauw(V1);
    dw2  = (P.winf(V2) -w2) / P.tauw(V2);
    dh1  = (P.hinf(V1) -h1) / P.tauh(V1);
    dm1  = (P.minf(V1) -m1) / P.taum(V1);
    dm2  = (P.minf(V2) -m2) / P.taum(V2);
    dh2  = (P.hinf(V2) -h2) / P.tauh(V2);
    dp   = (P.pinf(V1) - p) / P.taup(V1);
    da1  = (P.ainf(V1)- a1) / P.taua(V1);
    da2  = (P.ainf(V2)- a2) / P.taua(V2);
    %%%%% output %%%%%
    dx = [dV1 ; dV2 ; dw1 ; dh1; dw2; dm1;dm2; dh2; dp; da1; da2];
end 