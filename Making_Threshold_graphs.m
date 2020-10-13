%%%This file makes Threshold graphs for Sine, Ramp, Ramp2, and EPSGpair it
%%%calls BinarySearch.m


%same node for all runs (can be 3 or 5)
node=3;



%Combine_all = (stimType, number of points on graph (only change if you change graphing 
% bellow), node, factor away from Soma, max current value/voltage for
% EPSGpair)
%that classifies a spike)
[Thresholds1_ramp,Thresholds2_ramp]=BinarySearch('ramp', 15,node, 30, 15000);
[Thresholds1_sine,Thresholds2_sine]=BinarySearch('sine', 12,node, 30, 15000);
[Thresholds1_ramp2,Thresholds2_ramp2]=BinarySearch('ramp2', 15,node, 30, 15000);
[Thresholds1_EPSGpair,Thresholds2_EPSGpair]=BinarySearch('EPSGpair', 25,node,10, 150);



%%%Graphing begins...

%graph for ramp
Xvar_ramp=zeros(15,1);
for i=1:15
    Xvar_ramp(i)=5+(i/10);
end 
figure(1)
plot(Xvar_ramp, [Thresholds1_ramp,Thresholds2_ramp],'linewidth',2);
title(['Ramp Thresholds Compartment', num2str(node)])
xlabel('Stop Value')
ylabel('Threshold Voltage');
legend({'Multi-Compartment','Two-Compartment'},'fontsize',8,'box','off');

%graph for sine
Xvar_sine= zeros(12,1);
for i=1:12
    Xvar_sine(i)= i*100;
end 
figure(2)
plot(Xvar_sine, [Thresholds1_sine,Thresholds2_sine],'linewidth',2);
title(['Sine Thresholds Compartment', num2str(node)])
xlabel('Frequency')
ylabel('Threshold Voltage');
legend({'Multi-Compartment','Two-Compartment'},'fontsize',8,'box','off');


%graph for ramp2
Xvar_ramp2=zeros(15,1);
for i=1:15
    Xvar_ramp(i)=5+(i/10);
end 
figure(3)
plot(Xvar_ramp, [Thresholds1_ramp2,Thresholds2_ramp2],'linewidth',2);
title(['Ramp2 Thresholds Compartment', num2str(node)])
xlabel('Stop Value')
ylabel('Threshold Voltage');
legend({'Multi-Compartment','Two-Compartment'},'fontsize',8,'box','off');


%graph for EPSGpair
Xvar_EPSGpair=zeros(25,1);
for i=1:25
    Xvar_EPSGpair(i)= 5+((i-1)/25);
end 
figure(4)
plot(Xvar_EPSGpair, [Thresholds1_EPSGpair,Thresholds2_EPSGpair],'linewidth',2);
title(['EPSGpair Thresholds Compartment', num2str(node)])
xlabel('Time Difference')
ylabel('Threshold Voltage');
legend({'Multi-Compartment','Two-Compartment'},'fontsize',8,'box','off');



