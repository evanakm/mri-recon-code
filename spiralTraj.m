function [trajectory, kres] = spiralTraj(kReach,numSteps,numLoops,numLeaves,FOV)

    % kReach = extent, should be at least res/2

    kres = 2*numLoops;
    
    dTH = 2*pi*numLoops/numSteps;
%    x = 0:dTH:(2*pi*numLoops);
    x = 0:dTH:(2*pi*numLoops);
    x = x(2:end);
    nx = length(x);
    trajectory = zeros(nx*numLeaves,2);
    
    current = 1;
    for i = 1:numLeaves
        initPhase = 2*pi*i/numLeaves;
        drops = exp(1i*(x+initPhase))/FOV;
        trajectory(current:(current+nx-1),1) = x.*kReach/(2*pi*numLoops).*real(drops);
        trajectory(current:(current+nx-1),2) = x.*kReach/(2*pi*numLoops).*imag(drops);
        
        current = current+nx;

    end
%    scatter(trajectory(:,1),trajectory(:,2));

end