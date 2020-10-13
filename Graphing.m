function Graphing (graph, t1,y,t2,x,t3,z)
clf %Can be deleted when needed for runtime

if nargin>3
    V1 = x(:,1); V2 = x(:,2);
end
U1 = y(:,1); U2 = y(:,graph.node);

switch(graph.Type)
    case('Contour')
        [X,Y] = meshgrid(1:45,t1)
        surf(X,Y,y(:,1:45),'linestyle','none')
        title('Voltage Propogation')
        xlabel('Compartment Number (1-45)')
        ylabel('Time');
        zlabel('Voltage (mV)')
    
    case('ModelComparison')
        subplot(2,1,1)
             plot(t1,y(:,[ 1 graph.node]),'linewidth',2)
             title('multi-compartment')
             set(gca,'fontsize',12); 
             legend({'Soma',['node ',num2str(graph.node)]},'fontsize',8,'box','off'); 
             xlabel('Time (ms)')
             ylabel('Voltage (mV)');
             xlim([0 graph.tEnd])


        subplot(2,1,2)
            hold all; 
            plot(t2,[V1 V2],'linewidth',2)
            title('two-compartment')
            set(gca,'fontsize',12); 
            legend({'Soma','node'},'fontsize',8,'box','off'); 
            xlabel('Time (ms)')
            ylabel('Voltage (mV)');
            xlim([0 graph.tEnd])

        sgtitle('Comparison of 2 Compartment and Multi-Compartment Models')
    case('ModelComparison2')
    Z1 = z(:,1); Z2 = z(:,2);
    subplot(2,1,1)
        hold all; 
        plot(t2,[V1 V2],'linewidth',2)
        title('Two-Compartment (Node 3)')
        set(gca,'fontsize',12); 
        legend({'Soma','node'},'fontsize',8,'box','off'); 
        xlabel('Time (ms)')
        ylabel('Voltage (mV)');
        xlim([0 graph.tEnd])

    subplot(2,1,2)
        hold all; 
        plot(t3,[Z1 Z2],'linewidth',2)
        title('Two-Compartment (Node 5)')
        set(gca,'fontsize',12); 
        legend({'Soma','node'},'fontsize',8,'box','off'); 
        xlabel('Time (ms)')
        ylabel('Voltage (mV)');
        xlim([0 graph.tEnd])
        sgtitle('Comparison of 2 Compartment Models')
    case('NodeComparison')

        subplot(2,1,1)
            hold all
            plot(t1,U1,'linewidth',2)
            plot(t2,V1,'linewidth',2)
            title('Soma')
            set(gca,'fontsize',12); 
            legend({'Multi-Compartment','2 Compartment'},'fontsize',8,'box','off'); 
            xlabel('Time (ms)')
            ylabel('Voltage (mV)');
            xlim([0 graph.tEnd])

        subplot(2,1,2)
            hold all; 
            plot(t1,U2,'linewidth',2)
            plot(t2,V2,'linewidth',2)
            title(['node ',num2str(graph.node)])
            set(gca,'fontsize',12); 
            legend({'Multi-Compartment','2 Compartment'},'fontsize',8,'box','off');
            xlabel('Time (ms)')
            ylabel('Voltage (mV)');
            xlim([0 graph.tEnd])

         sgtitle('Comparison of Soma and graph.node')

    case('MultiCompartment')
        plot(t1,y(:,[ 1 graph.node]),'linewidth',2) %TODO w,z for Soma and graph.node compare like graph.nodecomparison
        title(['Multi-compartment model, node:', num2str(graph.node)])
        set(gca,'fontsize',12); 
        legend({'Soma',['node ',num2str(graph.node)]},'fontsize',8,'box','off'); 
        xlabel('Time (ms)')
        ylabel('Voltage (mV)');
        xlim([0 graph.tEnd])

    case('TwoCompartment')
        V1 = y(:,1); V2 = y(:,2);
        plot(t1,[V1 V2],'linewidth',2)
        title('two-compartment')
        set(gca,'fontsize',12); 
        legend({'Soma','node'},'fontsize',8,'box','off'); 
        xlabel('Time (ms)')
        ylabel('Voltage (mV)');
        xlim([0 graph.tEnd])
        
    case('Input')
        
                %%%% Setup Variables %%%%
        stop = graph.stop; I = graph.I; Syn = graph.Syn;
        Iext = zeros(size(t1)); stimType = graph.stimType ;inputNode = graph.inputNode;
                %%%% Caluclates Input %%%%
        
    switch(stimType); case('Synaptic');[Syn.tSyn, Syn.gSyn] = Synaptic(Syn);end
        
    for T=1:size(t1)

        switch(stimType)
        case('ramp')
            startRamp = graph.start;    % time at graph.start of ramp
            stopRamp = stop;     % slope of ramp [pA/ms]
            IrampMax  = I; % maximum current of ramp
            ramp= 1/(stopRamp-startRamp);
            tTop = startRamp+(1 / ramp);
            I0 = ((t1(T)>=startRamp)*(t1(T)<=tTop).*I.*(t1(T)-startRamp)*ramp);
            if t1(T)>graph.start && t1(T)<=graph.tEnd; Iext(T) =I0 ; end
        case('ramp2')
            startRamp = graph.start;    % time at graph.start of ramp
            stopRamp = stop;     % slope of ramp [pA/ms]
            IrampMax  = I; % maximum current of ramp
            ramp2= 1/(stopRamp-startRamp);
            tTop = startRamp+(1 / ramp2);
            Ia = ((t1(T)>=startRamp)*I*(t1(T)-startRamp)*ramp2)/1000;
            Ib = I/1000;
            I0 = min(Ia,Ib);
            Iext(T) =I0;
        case('step')
            if t1(T)>graph.start && t1(T)<=stop; Iext(T) = I; end
         case('sine')
            I0 = I(t1(T));
            if t1(T)>graph.start && t1(T)<= stop; Iext(T) = I0; end 
         case('Synaptic')
            V = y(:,inputNode);
            q = find(Syn.tSyn>=t1(T),2,'first');
            g = interp1q(Syn.tSyn,Syn.gSyn',t1(T));
            Iapp =  g * (V(T) - Syn.VsynE); 
            if t1(T)>graph.start && t1(T)<= stop; 
                Iext(T) = -Iapp; end
        case('EPSG')
            V = x(:,inputNode);
            tEPSG = t1(T)-graph.start;
            epsg = (1/0.21317) .* (tEPSG>0) .* (exp(-tEPSG/.18) - exp(-tEPSG/.1) ); % unitary epsg waveform
            Iapp = (I.*(0-V(T)).*(tEPSG>=0).*epsg) ;
            if t1(T)>graph.start && t1(T)<= stop; Iext(T) = Iapp; end
        case('EPSGpair')
            V = x(:,inputNode);
            tEPSG = t1(T)-graph.start;
            startEPSG1 = graph.start;   % graph.start of first EPSG event
            startEPSG2 = stop; % graph.start of second EPSG event
            td = startEPSG2-startEPSG1;
            epsg = @(t,G) (G/0.21317) .* ones(size(t)) .* (t>0) .* (exp(-t/.18) - exp(-t/.1) ); % unitary epsg waveform
            Iapp = (I.*(0-V(T).* ((tEPSG>=0).*epsg(tEPSG,1) + ...
                (tEPSG>=(td)).*epsg(tEPSG-td,1))));          
            Iext(T)=Iapp;
        end
        end

        plot(t1,Iext,'linewidth',2)
        title('Input')
        set(gca,'fontsize',12); 
        xlabel('Time (ms)')
        ylabel('Current (pA)');
        xlim([0 graph.tEnd])
                
    case('Multigraph')
        Z1 = z(:,1); Z2 = z(:,2);
        
    subplot(2,3,[1 2 3])
        plot(t1,y(:,[ 1 3 5 9 13 21 29 37 45]),'linewidth',2)
        title('multi-compartment')
        set(gca,'fontsize',12); 
        legend({'Soma','node 3','node 5','node 9','node 13','node 21'...
            ,'node 29','node 37', 'node 45'},'fontsize',8,'box','off'); 
        xlabel('Time (ms)')
        ylabel('Voltage (mV)');
        xlim([0 graph.tEnd])

    subplot(2,3,4)
        hold all; 
        plot(t2,[V1 V2],'linewidth',2)
        title('Two-Compartment (Node 3)')
        set(gca,'fontsize',12); 
        legend({'Soma','node'},'fontsize',8,'box','off'); 
        xlabel('Time (ms)')
        ylabel('Voltage (mV)');
        xlim([0 graph.tEnd])

    subplot(2,3,5)
        hold all; 
        plot(t3,[Z1 Z2],'linewidth',2)
        title('Two-Compartment (Node 5)')
        set(gca,'fontsize',12); 
        legend({'Soma','node'},'fontsize',8,'box','off'); 
        xlabel('Time (ms)')
        ylabel('Voltage (mV)');
        xlim([0 graph.tEnd])
        %%%% Setup Variables %%%%
        stop = graph.stop; I = graph.I; Syn = graph.Syn;
        Iext = zeros(size(t1)); stimType = graph.stimType ;inputNode = graph.inputNode;
        
        %%%% Caluclates Input %%%%
        
    switch(stimType); case('Synaptic');[Syn.tSyn, Syn.gSyn] = Synaptic(Syn);end
        
    for T=1:size(t1)

        switch(stimType)
        case('ramp')
            startRamp = graph.start;    % time at graph.start of ramp
            stopRamp = stop;     % slope of ramp [pA/ms]
            IrampMax  = I; % maximum current of ramp
            ramp= 1/(stopRamp-startRamp);
            tTop = startRamp+(1 / ramp);
            I0 = ((t1(T)>=startRamp)*(t1(T)<=tTop).*I.*(t1(T)-startRamp)*ramp);
            if t1(T)>graph.start && t1(T)<=graph.tEnd; Iext(T) =I0 ; end
        case('ramp2')
            startRamp = graph.start;    % time at graph.start of ramp
            stopRamp = stop;     % slope of ramp [pA/ms]
            IrampMax  = I; % maximum current of ramp
            ramp2= 1/(stopRamp-startRamp);
            tTop = startRamp+(1 / ramp2);
            Ia = ((t1(T)>=startRamp)*I*(t1(T)-startRamp)*ramp2)/1000;
            Ib = I/1000;
            I0 = min(Ia,Ib);
            Iext(T) =I0;
        case('step')
            if t1(T)>graph.start && t1(T)<=stop; Iext(T) = I; end
         case('sine')
            I0 = I(t1(T));
            if t1(T)>graph.start && t1(T)<= stop; Iext(T) = I0; end 
         case('Synaptic')
            V = y(:,inputNode);
            q = find(Syn.tSyn>=t1(T),2,'first');
            g = interp1q(Syn.tSyn,Syn.gSyn',t1(T));
            Iapp =  g * (V(T) - Syn.VsynE); 
            if t1(T)>graph.start && t1(T)<= stop; 
                Iext(T) = -Iapp; end
        case('EPSG')
            V = x(:,inputNode);
            tEPSG = t1(T)-graph.start;
            epsg = (1/0.21317) .* (tEPSG>0) .* (exp(-tEPSG/.18) - exp(-tEPSG/.1) ); % unitary epsg waveform
            Iapp = (I.*(0-V(T)).*(tEPSG>=0).*epsg) ;
            if t1(T)>graph.start && t1(T)<= stop; Iext(T) = Iapp; end
        case('EPSGpair')
            V = x(:,inputNode);
            tEPSG = t1(T)-graph.start;
            startEPSG1 = graph.start;   % graph.start of first EPSG event
            startEPSG2 = stop; % graph.start of second EPSG event
            td = startEPSG2-startEPSG1;
            epsg = @(t,G) (G/0.21317) .* ones(size(t)) .* (t>0) .* (exp(-t/.18) - exp(-t/.1) ); % unitary epsg waveform
            Iapp = (I.*(0-V(T).* ((tEPSG>=0).*epsg(tEPSG,1) + ...
                (tEPSG>=(td)).*epsg(tEPSG-td,1))));          
            Iext(T)=Iapp;
        end
        end
      subplot(2,3,6)
        hold all;   
        plot(t1,Iext,'linewidth',2)
        title('Input')
        set(gca,'fontsize',12); 
        xlabel('Time (ms)')
        ylabel('Voltage (mV)');
        xlim([0 graph.tEnd])
        sgtitle('Comparison of 2 Compartment and Multi-Compartment Models with Input')

end
    
%%measures z1 value
if graph.z1 == true
    figure(2);
    clf;
    z1 = y(:,226);
    fitline = @(x,P) P(1) + P(2)*x;
    x = [0:graph.graph.tEnd];
    [fitP, ~, exitFlag]= fminsearch( @(P) norm( z1 -fitline(x,P)), [2.1 0]);
    plot(x,fitline(x,fitP),'k','linewidth',3)
    title('z1 vs time')
    disp(['z1=', num2str(fitP(1))]);
end

%%measures a value
if graph.a == true
    figure(3);
    clf;
    a = y(:,271);
    fitline = @(x,P) P(1) + P(2)*x;
    x = [0:graph.graph.tEnd];
    fitP = fminsearch( @(P) norm( a -fitline(x,P)), [2.1 0]);
    plot(x,fitline(x,fitP),'k','linewidth',3)
    title('a vs time')
    disp(['a=', num2str(fitP(1))]);
end

end