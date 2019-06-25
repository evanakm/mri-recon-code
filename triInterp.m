function [kPTS, kernK, kerntraj, gridsize] = triInterp(trajectory,res,FOV,kspace)

counter = size(trajectory,1);

kPTS = zeros(counter,1);
%kFOVx = abs(min(trajectory(:,1)));
%kFOVy = abs(min(trajectory(:,2)));

basePTS = zeros(counter,1);
%indices = zeros(counter,2); %Not necessary. Helpful for debugging.

centreOfFFT = res/2;


for i = 1:counter
    idxKX = trajectory(i,1)*FOV + centreOfFFT;
    idxKY = trajectory(i,2)*FOV + centreOfFFT;
    
%    indices(i,:) = [ idxKX idxKY]; %Not necessary. Only used for debugging.
    fracX = idxKX - floor(idxKX);
    fracY = idxKY - floor(idxKY);
    
    if(idxKX < -0.5 || idxKX >= res+0.5 || idxKY < -0.5 || idxKY >= res+0.5)
        kPTS(i) = 0; basePTS(i)=0;
    else
        
        % Convolving with a triangle
        if(fracX >= 0.5)
            lambdaX = 2 * fracX - 1; % reduces from 1 - (1-fracX)/.5, but this is not clear at first glance.
            xPT = ceil(idxKX);
        else
            lambdaX = 1 - 2 * fracX;
            xPT = floor(idxKX);
        end
        
        if(fracY >= 0.5)
            lambdaY = 2 * fracY - 1; % reduces from 1 - (1-fracX)/.5, but this is not clear at first glance.
            yPT = ceil(idxKY);
        else
            lambdaY = 1 - 2 * fracY;
            yPT = floor(idxKY);
        end
        
        if(xPT == 0)
            xPT = res;
        end
        
        if(yPT == 0)
            yPT = res;
        end
            
        kPTS(i) = lambdaX*lambdaY*kspace(xPT,yPT);
        basePTS(i) = lambdaX*lambdaY; % Convolving a grid of 1's.
    
    
    end
    
end
    
    % Inverting the convolutional kernel.
    % The grid goes from -(kres/2)*2pi/res/FOV to (kres/2)*2pi/res/FOV. The
    % "normal" grid has kres = res, and kres steps.
    DX = pi/res/FOV;
    [tx,ty] = meshgrid(linspace(-2*DX,2*DX,res*3),linspace(-2*DX,2*DX,res*3));
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

