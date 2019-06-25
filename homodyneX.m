function image = homodyneX(kspace,midstrip,filt)

    % ------- Dropped ------ | 2* midstrip | --------- Kept ---------
    
    [sz1,sz2] = size(kspace);
    
    if (midstrip > sz1/2) 
        error('midstrip too big\n');
    end
    
    if (nargin == 1)
        error('not enough arguments\n');
    elseif(nargin == 2)
        filt = 'step';
    end
    
    midpoint = floor(sz1/2);
    
    low = 1:(midpoint-midstrip);
    middle = (midpoint-midstrip+1):(midpoint+midstrip);
    high = (midpoint+midstrip+1):sz1;
    
    mask = zeros(sz1,sz2);
    lowpass = kspace;
    lowpass(low,:) = 0;
    lowpass(high,:) = 0;

    if (strcmp(filt,'step'))
        mask(middle,:) = 1;
        mask(high,:) = 2;
    elseif (strcmp(filt,'ramp'))
        [mesh1,mesh2] = meshgrid(ones(sz1,1),linspace(0,2,length(middle)));
        mask(middle,:) = mesh2;
        mask(high,:) = 2;
    else
        error('unknown filter\n');
    end
    
    ksp = mask .* kspace;
    
    refArr = ifft2(lowpass); %Lowpass
    Arr = ifft2(ksp);
    
    argu = angle(refArr);
    demod = Arr .* exp(-1i*argu);
    
    image = real(demod);
    % Let's invert
    
    
end

