function image = griddingAlgo(kspace,options)
    
    res     = options.res;
    interp  = options.interp;
    traj    = options.traj;
    kReach  = options.kReach;
    FOV     = 1; %The code base I'm using requires this parameter, but the following requires it to be kept at 1.
    tau     = options.tau;
    nearest = options.nearest;
    over    = options.over;
    

    switch traj
        case 'cartesian'
            fsamp = options.par1;
            R = options.par2;
            [traj, kres] = cartesianTraj(kReach,fsamp,R,FOV);
        case 'polar'
            rSteps = options.par1;
            thSteps = options.par2;
            [traj, kres] = polarTraj(kReach,rSteps,thSteps,FOV,true);
        case 'spiral'
            numSteps = options.par1;
            numLoops = options.par2;
            numLeaves = options.par3;
            [traj, kres] = spiralTraj(kReach,numSteps,numLoops,numLeaves,FOV);
        case 'custom'
            traj = options.par1;
            kres = res;
        otherwise
            error('unknown trajectory');
    end
    


    

    switch interp
        case 'linear'
            image = griddingTriangular(kspace,traj,res,kres,nearest,FOV);
        case 'gaussian'
            bottom = -res/2; top = (res/2)-1;
            [kXgrid,kYgrid] = meshgrid((bottom:top)/FOV,(bottom:top)/FOV);
            image = griddingGaussian(kspace,traj,res,kres,kXgrid,kYgrid,tau,nearest,FOV);
            image = image';
        otherwise
            error('unknown interpolation method');
    end
    
end
