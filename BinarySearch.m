
%%%passes back Thresholds1= Thresholds for 45-CPT model and Thresholds2 =
%%%Thresholds for 2-CPT model
%%% stimType, step_size, node, and factor are defined in the
%%% Making_Threshold_graphs.m


function [Thresholds1, Thresholds2]=BinarySearch(stimType, step_size, node, factor, max)  
type='active-full'; %Options:passive and active-KLT and active-sodium
node=node; %node in question (other than Soma)
inputNode = 1; % Any node that isn't 1 will result in an input into compartment 2
stimType=stimType; %current type
tEnd=20; % End of Data gathering
v0=-68; %base voltage (homeostasis)
start = 5;   % starting value of the input current
zoom=1; %how close we want the bonary search to look, don't do more than 100
step_size=step_size; %how many times we go through the binary search
max=max; %the max threshold we care about
factor=factor; %how far the CPT from Soma during a spike

%%%Threshold for MsoAxon
Thresholds1=zeros(step_size,1);

for i = 1:step_size
    
    first_t1= 0;
    location1=max;
    previous1=0;
    distance1=max;
    
    %standard binary search algorithm follows....
    while distance1>zoom  && location1<=max

        I = location1; %Input Current or Max current
        switch(stimType)

            case('step')
                stop  =  15;  % time at end of step

            case('ramp')
                stop = 5+(i/10);     % slope of ramp [pA/ms]

            case('ramp2')  %ramp doesn't go back down
                stop = 5+(i/10);

            case('sine')
                Syn.f = i * 100;      %frequency of wave        
                stop = start + 1000 /(2*Syn.f);   % end val don't change

            case('EPSG')                        
                stop = 10; %Ignore, unused

            case('EPSGpair')
                
                start2 = start+((i-1)/25); % start of second EPSG event 
                
                stop = start2; %ignore

        end
        
        
        switch(stimType)
            case('Synaptic')
                Syn.tEnd = tEnd; %simulation length
                Syn.freq = 350; %input frequency (Hz)
                Syn.gE   = 1.3 * 1e7; % epsg conductance (non-meaningful units at the moment)
                Syn.VsynE = 0;   % reversal potential for excitation (mV). not given in paper?
                Syn.randomIn = 13986; %Random Input Syn

                I = 0;
                stop = 15;
                Syn.Exists = true;
            otherwise
                Syn.Exists = false;
        end

        [t1,y] = msoAxon(stimType,start,stop,I,node,type,tEnd,v0,inputNode, Syn);
        
        %%%Spike Detection for msoAxon
        spike_count1=0;
        factor1=factor; %amount of Node above Soma
        re_set1=0; %make sure we don't count a spike again
        U1 = y(:,1); 
        U2 = y(:,node);

        for l=1:length(t1)    
            if (U2(l)-U1(l)>factor1) && re_set1==0        
                spike_count1=spike_count1+1; %if the CPT 2 is in spiking range, a spike is recorded
                re_set1=1;             
            end 
            if U2(l)<U1(l)
                re_set1=0; %another spike can only be recorded once CPT 2 goes back to resting, below the Soma      
            end
        end

        dummy_location1=location1; %so we can write over location1       
        first_t1=location1; %location between where we know a spike occured and one didn't
        distance1=abs((location1-previous1)/2); %distance between the current location and the previously tested
        

        if spike_count1~=0        
            location1=location1-distance1; %if there was a spike we must go closer to 0 to find a nonspiking region     
        end 

        if spike_count1==0        
            location1=location1+distance1; %if there was a spike we must go closer to max to find a spiking region

        end 
        previous1=dummy_location1; %our location1 has been written over, this is where dummlocation1 is employed
    
    end 
    
   switch( stimType)
        case('EPSGpair')
            Thresholds1(i)=round(first_t1,0); %we round to the 1st digit with EPSGpair since our max is 150
        otherwise 
            Thresholds1(i)=round(first_t1,-1); %we round to the 10s digit since our max is 15000
    end 
end 


%%%Threshold for EasyRun

Thresholds2=zeros(step_size,1);

for i=1:step_size%step-1, ramp-15, ramp2-15, sine-12,EPSG-1, EPSGpair-11
    
    first_t2= 0;
    location2=max;
    previous2=0;
    distance2=max;
    
    %standard binary search algorithm follows...
    while distance2>zoom && location1<=max
        I = location2; %Input Current or Max current 

        switch(stimType)

            case('step')
                stop  =  15;  % time at end of step

            case('ramp')
                stop = 5+(i/10);     % slope of ramp [pA/ms]

            case('ramp2')  %ramp doesn't go back down
                 stop = 5+(i/10);

            case('sine')
                Syn.f = i * 100;      %frequency of wave    
                stop = start + 1000 /(2*Syn.f);   % end val don't change

            case('EPSG')
                stop = 10; %Ignore, unused

            case('EPSGpair')
                %I = 70;%I/400; % amplitude of EPSG event [26.7 is a "unitary EPSG" in the manuscript]
                start2 = start+((i-1)/25); % start of second EPSG event        
                stop = start2; %ignore

        end
        
        switch(stimType)
            case('Synaptic')
                Syn.tEnd = tEnd; %simulation length
                Syn.freq = 350; %input frequency (Hz)
                Syn.gE   = 1.3 * 1e7; % epsg conductance (non-meaningful units at the moment)
                Syn.VsynE = 0;   % reversal potential for excitation (mV). not given in paper?
                Syn.randomIn = 13986; %Random Input Syn

                I = 0;
                stop = 15;
                Syn.Exists = true;
            otherwise
                Syn.Exists = false;
        end


        [t2,x] = TwoCpt(stimType,start,stop,I,node,type,tEnd,v0,inputNode, Syn);

        
        V1 = x(:,1); 
        V2 = x(:,2);

        %%%Spike Detection for msoAxon
        spike_count2=0;
        factor2=factor; %amount of Node above Soma
        re_set2=0; %make sure we don't count a spike again   
        for l=1:length(t2)    
            if (V2(l)-V1(l)>factor2) && re_set2==0        
                spike_count2=spike_count2+1;%if the CPT 2 is in spiking range, a spike is recorded
                re_set2=1;             
            end 
            if V2(l)<V1(l)
                re_set2=0;  %another spike can only be recorded once CPT 2 goes back to resting, below the Soma      
            end
        end 

        dummy_location2=location2; %so we can write over location2
        distance2=abs((location2-previous2)/2); %distance between the current location and the previously tested       
        first_t2=location2; %location between where we know a spike occured and one didn't

        if spike_count2~=0        
            location2=location2-distance2;      
        end 

        if spike_count2==0        
            location2=location2+distance2;

        end 
        previous2=dummy_location2; 
    end
    
    
    switch( stimType)
        case('EPSGpair')
            Thresholds2(i)=round(first_t2,0);
        otherwise 
            Thresholds2(i)=round(first_t2,-1);
    end 
    
end