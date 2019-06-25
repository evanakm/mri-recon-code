function [kPTS, kernK, kerntraj, gridsize] = fullInterp(trajectory,res,FOV,kspace)

counter = size(trajectory,1);

kPTS = zeros(counter,1);
%kFOVx = abs(min(trajectory(:,1)));
%kFOVy = abs(min(trajectory(:,2)));

basePTS = zeros(counter,1);
indices = zeros(counter,2); %Not necessary. Helpful for debugging.

centreOfFFT = res/2;


for i = 1:counter
    idxKX = trajectory(i,1)*FOV + centreOfFFT; %Re-indexing it for squareinterp
    idxKY = trajectory(i,2)*FOV + centreOfFFT;
    
    indices(i,:) = [ idxKX idxKY];
    
    if(idxKX < 0 || idxKX >= res+1 || idxKY < 0 || idxKY >= res+1)
        kPTS(i) = 0; basePTS(i)=0;        
    else
   
        if(idxKX >= 0 && idxKX < 1)
            idxKX = idxKX + res;
        end
        if(idxKY >= 0 && idxKY < 1)
            idxKY = idxKY + res;
        end
    
        kPTS(i) = squareinterp(idxKX,idxKY,kspace);
%        basePTS(i) = squareinterp(idxKX,idxKY,kbase);
    end
    

end
    
    % Inverting the convolutional kernel.
    % The grid goes from -(kres/2)*2pi/res/FOV to (kres/2)*2pi/res/FOV. The
    % "normal" grid has kres = res, and kres steps.
    DX = 2*pi/res/FOV;
    [tx,ty] = meshgrid(linspace(-2*DX,2*DX,res*3),linspace(-4*DX,4*DX,res*3));
    kernel = @(x,a) max(1 - abs(x/a),0);
    GRID = kernel(tx,DX).*kernel(ty,DX);
    
    kernsz = numel(GRID);
    kerntraj = zeros(kernsz,2);
    kernK = zeros(kernsz,1);
    
    for i=1:kernsz
        kerntraj(i,1) = tx(i);
        kerntraj(i,2) = ty(i);
        kernK(i) = GRID(i);
    end

    gridsize = size(GRID,1);
    
    
end

