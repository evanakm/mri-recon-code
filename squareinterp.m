function inter = squareinterp(x,y,kspace)
    xfl = floor(x);
    yfl = floor(y);
    lambda = x - xfl; %decimal part of x
    mu = y - yfl; %decimal part of y

    % edge case handling
    
    xmax = size(kspace,1);
    ymax = size(kspace,2);
    
    if(xfl == xmax && yfl == ymax)
        step1 = kspace(xfl,yfl) + lambda*(kspace(1,yfl) - kspace(xfl,yfl));
        step2 = kspace(xfl,1) + lambda*(kspace(1,1) - kspace(xfl,1));
        inter = step1 + mu*(step2-step1);
    elseif(xfl==xmax)
        step1 = kspace(xfl,yfl) + lambda*(kspace(1,yfl) - kspace(xfl,yfl));
        step2 = kspace(xfl,yfl+1) + lambda*(kspace(1,yfl+1) - kspace(xfl,yfl+1));
        inter = step1 + mu*(step2-step1);
    elseif(yfl==ymax)
        step1 = kspace(xfl,yfl) + lambda*(kspace(xfl+1,yfl) - kspace(xfl,yfl));
        step2 = kspace(xfl,1) + lambda*(kspace(xfl+1,1) - kspace(xfl,1));
        inter = step1 + mu*(step2-step1);
    else
        step1 = kspace(xfl,yfl) + lambda*(kspace(xfl+1,yfl) - kspace(xfl,yfl));
        step2 = kspace(xfl,yfl+1) + lambda*(kspace(xfl+1,yfl+1) - kspace(xfl,yfl+1));
        inter = step1 + mu*(step2-step1);
    end
    
end
