function [convol, kernK, kerntraj, gridsize] = gaussianInterp(trajectory,kspace,kXgrid,kYgrid,tau,res,FOV)

kpts = size(trajectory,1);

gaussian = @(x,y,sigma) exp(-1*(x.^2+y.^2)/(2*sigma^2));
convol = zeros(kpts,1);

for i=1:kpts
    gaussGrid = gaussian(kXgrid-trajectory(i,1),kYgrid-trajectory(i,2),tau);
    convol(i) = sum(sum(gaussGrid.*kspace));
end

    % Inverting the convolutional kernel.
    % The grid goes from -(kres/2)*2pi/res/FOV to (kres/2)*2pi/res/FOV. The
    % "normal" grid has kres = res, and kres steps.
    realTau = 2*pi*tau/res;
    
    [tx,ty] = meshgrid(linspace(-3*realTau,3*realTau,res*3),linspace(-3*realTau,3*realTau,res*3));

    GRID = gaussian(tx,ty,realTau);
    
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