function Spikes = Spiking(x,factor, Model)
    
    switch(Model)
        case{'Two', 'two'}
        V1 = x(:,1);
        Diff = zeros(size(x,1), 2);
        for i = 1:2
            Diff(:,i) = x(:,i) -V1;
        end
        
    otherwise
        V1 = x(:,1);
        
        Diff = zeros(size(x,1), 45);
        for i = 1:45
            Diff(:,i) = x(:,i) -V1;
        end
    end

    Spikes = zeros(size(Diff));
    reset  = 0;
    
    for j = 1:size(Spikes,2)
    for i = 1:size(Spikes,1)
        if (Diff(i,j) > factor) && (reset == 0)
            Spikes(i,j) = 1;
            reset = 1;
        end
        if (Diff(i,j) <= 0)
            reset = 0;
        end
    end
    end   
end