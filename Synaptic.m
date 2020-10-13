function [tSyn,gSyn] = Synaptic(Syn)

    % outputs are time vector and conductance vector for synaptic input
    
    % white noise
    dt = 0.01;
    tSyn = 0:dt:Syn.tEnd;
    nt = length(tSyn);
    rng(Syn.randomIn);%seeds random number
    n = randn(nt,1);

    % center frequency
    wc = Syn.freq*2*pi/1000; 

    % convolve with gammatone filter
    gw = 24.7 * (4.37 * wc/(2*pi) + 1);
    gamFilter = (tSyn/1000).^4 .* exp(-2*pi*tSyn*gw/1000) .* cos(tSyn*wc);
    gamConv = conv(gamFilter,n); 
    
    % rectify
    rectGam = max( gamConv(1:nt), 0 );
    
    % convolve with epsg waveform
    % 17.57 normalizes so max is G, 37 is unitary epsg amp , as in Lehnert et al
    epsg = Syn.gE*37*17.57*( (1 - exp(-tSyn/1.)).^1.3  .* exp(-tSyn/0.27) ) ; 
    epsgConv = conv(rectGam,epsg);
    gSyn = epsgConv(1:nt);    
end