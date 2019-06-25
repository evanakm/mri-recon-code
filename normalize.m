    function [norm,maxim] = normalize(x)
        
        maxim = max(abs(x(:)));
        norm = x/maxim;
    end
