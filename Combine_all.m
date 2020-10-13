clear all;

%{

**stimType: The type of stimulus being fed into our model
options: (step,ramp, ramp2, sine, epsg, epsg2, synaptic)

**node: The node we want to recieve dat about, while works for all
nodes, nodes 3 and 5, which are the first AIS segment and first Internode.

**inputNode: The Node in which current is introduced, for Two CPt if it is
larger than 1 it will be introduced into the second compartment a valid
change could be made so that it is the larger than 2. The range is 1 to 45.
There is less confluence of models after the first 1 or 2 nodes.

**stimType: The type of 
%}

type='active-full'; %Options:passive and active-KLT and active-sodium
node=3; %node in question (other than Soma)
inputNode = 1; % Any node that isn't 1 will result in an input into compartment 2
stimType='EPSGpair'; % Problem including Synaptic pair to tEnd not sure why?
tEnd=20; % End of Data gathering
v0=-68;
graph.Type  = 'ModelComparison'; %Options: ModelComparison and NodeComparison
graph.z1 = false; %Measures the z gating variable
graph.a  = false; %Measures the z gating variable
factor=30; %30 for all but 10 for EPSG, amount above Soma that counts as a spike

start = 5;    % time at start of input
I = 8500; %Input Current or Max current

switch(stimType)
    case('step')
        stop  =  10;  % time at end of step
    case('ramp')
        stop = 5;     % slope of ramp [pA/ms]
    case('ramp2')%ramp followed by never step till end
        stop = 5.5;
    case('sine')
        Syn.f = 200;      %frequency of wave        
        stop = start + 1000 /(2*Syn.f);   % end val don't change
    case('EPSG')
        I = 108;%I/400; % amplitude of EPSG event [26.7 is a "unitary EPSG" in the manuscript]        
        stop = 10; %Ignore, unused
    case('EPSGpair')
        I = 40;
        start2 = 5; % start of second EPSG event        
        stop = start2; %ignore
end
switch(stimType)
    case{'Synaptic', 'SynapticPair'}
        Syn.tEnd = tEnd; %simulation length
        Syn.freq = 350; %input frequency (Hz)
        Syn.gE   = 1.3 * 1e7; % epsg conductance (non-meaningful units at the moment)
        Syn.VsynE = 0;   % reversal potential for excitation (mV). not given in paper?
        Syn.randomIn = 1396; %Random Input Syn
        Syn.Diff = 1; %time difference in synaptic pair, must be smaller than difference in stop and tEnd
        
        I = 0;
        stop = 15;
        Syn.Exists = true;
    otherwise
        Syn.Exists = false;
end

%%Runs Models
[t1,y]  = msoAxon(stimType,start,stop,I,node,type,tEnd,v0,inputNode, Syn);
[t2,x]  = TwoCpt(stimType,start,stop,I,node,type,tEnd,v0,inputNode, Syn);
[t3,z]  = TwoCpt(stimType,start,stop,I,node,type,tEnd,v0,inputNode, Syn);


%%Spike Detection TwoCpt
Spike1 = Spiking(x,factor, 'Two'); %Spiking in TwoCpt
spikeCountTwo = 0;
for i=2:length(Spike1)    
    if (Spike1(i,1) > 0) && (Spike1(i-1,1) <= 0)
        spikeCountTwo=spikeCountTwo+1;
    end
end

%%%Spike Detection MultCpt
Spike2 = Spiking(y,factor, 'Mult'); %Spiking in MultCpt
spikeCountMult=zeros(1,45);
for j = 1:45
    SpikeNode  = Spike2(:,j);
    for i=2:length(SpikeNode)
        if (SpikeNode(i) > 0) && (SpikeNode(i-1) <= 0)
            spikeCountMult(1,j)=spikeCountMult(1,j) +1;
        end
    end
end

%%Graphs models
graph.node = node; graph.tEnd = tEnd; graph.stimType = stimType ; graph.inputNode = inputNode;
graph.start = start; graph.stop = stop; graph.I = I; graph.Syn = Syn; 

Graphing (graph,t1,y,t2,x,t3,z);