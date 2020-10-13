Vrest=repmat(-68,45,1);

%%%H
gh = [0.02 ; 0.02 ; 0.02 ; zeros(42,1)];
%ainf = 1./ (1+exp(0.1*(V+80.4)));
%Ih = gh .*  a .* (V - Vh);
%therefore h is...
h=gh.*(1./ (1+exp(0.1*(Vrest+80.4))));
 
%%%NA
gNa = [0.2 ; 4 ; 4 ; repmat([0 ; 4], 21,1)];
%minf = 1./(1 + exp((V+46)/(-11)));
%hinf = 1 ./  (1 + exp((V+62.5)/7.77));
%INa  = gNa .* m.^4 .* (0.993*h +0.007) .* (V - VNa); 
%therefore NA is...
Na=gNa.*((1./(1 + exp((Vrest+46)/(-11)))).^4).*(((1./(1 + exp((Vrest+62.5)/7.77)))*.993)+.007);
 
%%%KHT
gKHT = zeros(45,1); gKHT(1) = 0.1;
%pinf = 1./ ( 1 + exp(-(V+23)/6));
%IKHT =gKHT .* p .* (V - VK);
%therefore KHT is...
KHT=gKHT.*(1./ ( 1 + exp(-(Vrest+23)/6)));

%%%KLT
gKLT = [1.55 ; 1.55 ; 1.55 ; repmat([0 ; 1.55],21,1)];
%winf = 1 ./ (1 + exp((V+57.34)/(-11.7)));
%zinf = (1-0.27)./(1 + exp( (V+67) / 6.16) ) + 0.27;
%IKLT = gKLT .* w.^4 .* z .* (V - VK);
%therefore KLT is...
KLT= gKLT.* ((1 ./ (1 + exp((Vrest+57.34)/(-11.7)))).^4).*((1-0.27)./(1 + exp( (Vrest+67) / 6.16) ) + 0.27);

%%%glk
glk = [0.0005 ; 0.0005 ; 0.0005 ; repmat([0.0002 ; 0.05],21,1)];
%Ilk=glk

%%% Absolute Total Conductance
Total=KLT+KHT+Na+h+glk;

%%% Calcutes KLT Frac
KLT_frac=zeros(45,1);
for i=1:45
    if Total(i)~=0
        KLT_frac(i)=KLT(i)/Total(i);
    end
end

%%% Calculates H Frac
H_frac=zeros(45,1);
for i=1:45
    if Total(i)~=0
        H_frac(i)=h(i)/Total(i);
    end
end

%%% Calcutes Na Frac
Na_frac=zeros(45,1);
for i=1:45
    if Total(i)~=0
        Na_frac(i)=Na(i)/Total(i);
    end
end

save('Fractions.mat', 'KLT_frac', 'H_frac','Na_frac', '-v7.3', '-nocompression')

%Areas Surface Areas Calculations

% membrane surface area [um2]
% using equation for surface area of conical frustrum for tapered AIS, using inner are for myelinated internodes
SA = [8750 ; pi*(1.64/2 + 0.66/2)*sqrt( (1.64/2-1.66/3)^2+ 10^2) ; pi*10*.66 ; repmat([pi*100*0.66 ; pi*1*0.66],21,1)]; 
SAcm = SA*1e-8; %[cm^2]
areaRatio = SAcm ./ SAcm(1); %% add as a constant elsewhere so isn't Recalculated Eachtime
    
% length of each compartment [um]. 
% for soma: instead of sphere, use cylinder with length same as diameter of sphere
L = [sqrt(8750/pi) ; 10 ; 10 ; repmat([100 ; 1], 21, 1) ]; 
Lcm = L*1e-4; % [cm]
    
% cross sectional area [um2]
% for soma: set radius set to get same volume as sphere.
% for tapered compartment, get rid of taper. say its a cylinder with radius set to keep same surface area
Rcyl = SA./(2*pi*L); % compartment radius, assuming all compartments are cylinders, using A = pi*(2r)*L 
XA = pi*Rcyl.^2;
XAcm = XA*1e-8; %[cm^2] 

save('Area.mat', 'SA', 'SAcm','areaRatio', 'L', 'Lcm', 'Rcyl', 'XA', 'XAcm', '-v7.3', '-nocompression')