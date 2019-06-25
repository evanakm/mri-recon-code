function [trajectory, kres] = polarTraj(kReach,rSteps,thSteps,FOV,skipzero)

% kReach    how far it reaches in k-space, multiply by 2pi/res/FOV to get
%           appropriate units.
% rSteps    how many steps taken in the r direction. rSteps = res and
%           kReach = res/2 give the Nyquist.


eps = 0.01;

%kReach = res/2;
%rSteps = res;

dr = 2*kReach/rSteps;
RR = 0:dr:kReach;
centreOfkSpace = length(RR);

RR = [ -RR(end:-1:2) RR ];
RR = RR(1:(end-1));

kres = length(RR);

dth = pi/(thSteps+1);

TH = 0:dth:(pi-eps); %Makes sure we don't quite get to pi.


hitZeroYet = false; %Makes sure we only get one point at the origin, instead of many repeated points
counter = 0;

trajectory = zeros(length(RR)*length(TH),2);
kPTS = zeros(size(trajectory,1),1);

for th = TH
    for r = RR

        spot = r * exp(1i*th);
            
        if(abs(r)>0.001 || ~skipzero)
            counter = counter + 1;
            trajectory(counter,1) = real(spot)/FOV;
            trajectory(counter,2) = imag(spot)/FOV;
        elseif(~hitZeroYet && skipzero)
            counter = counter + 1;
            trajectory(counter,1) = real(spot)/FOV;
            trajectory(counter,2) = imag(spot)/FOV;
            hitZeroYet = true;
        end
        
        
    end
end

trajectory = trajectory(1:counter,:);

