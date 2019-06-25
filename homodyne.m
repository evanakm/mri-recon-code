function homodyne(kspace,midstrip,rowsORcols,filt)

    % ------- Dropped ------ | 2* midstrip | --------- Kept ---------

    dims = size(kspace);
    
    if (any(midstrip > dims/2)) 
        error('midstrip too big\n');
    end
    
    if (nargin < 2)
        error('not enough arguments\n');
    elseif(nargin == 3)
        filt = 'zero-one';
    end
    
    low = 1:
    indices = midpoint-midstrip+1:midpoint+midstrip;
    
    if (strcmp(filt,'zero-one'))
        mask = zeros(dims);
        filt = ones(2*midstrip,size(kspace,2));
    elseif (strcmp(filt,'ramp'))
        midpoint = floor(filt/2);
        [mesh1 mesh2] = meshgrid(,
        
    else
        error('unknown filter\n');
    end
    
    
    kspace
    
end

