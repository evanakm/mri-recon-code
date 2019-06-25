function image = halfKspace(kspace,direction)

    switch direction
        case 'kx'
            res = size(kspace,2);
            mid = res/2;
   
            Ksp = zeros(size(kspace));
            Ksp(:,1:(mid+1)) = kspace(:,1:(mid+1));
            Ksp(:,(mid+2):res) = conj(kspace(res:-1:1,mid:-1:2));
            
            %need to do one final realignment
            
            lastRow = Ksp(res,(mid+2):res);

            Ksp(2:res,(mid+2):res) = Ksp(1:(res-1),(mid+2):res);
            Ksp(1,(mid+2):res) = lastRow;
        case 'ky'
            res = size(kspace,1);
            mid = res/2;
   
            Ksp = zeros(size(kspace));
            Ksp(1:(mid+1),:) = kspace(1:(mid+1),:);
            Ksp((mid+2):res,:) = conj(kspace(mid:-1:2,res:-1:1));

            %need to do one final realignment
            
            lastRow = Ksp((mid+2):res,res);

            Ksp((mid+2):res,2:res) = Ksp((mid+2):res,1:(res-1));
            Ksp((mid+2):res,1) = lastRow;
            
        otherwise
            error('must be kx or ky');
    end
        
    image = ifft2(Ksp);
    
end