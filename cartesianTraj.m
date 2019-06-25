function [trajectory, kres] = cartesianTraj(kspan,fsamp,R,FOV)

% kspan    how far it reaches in k-space, multiply by 2pi/res/FOV to get
%           appropriate units.
% fsamp     steps taken in the kx direction.
% R         steps taken in the ky direction.


trajectory = zeros(floor(kspan/fsamp)*floor(kspan/R),2);

counter = 0;
for i=0:R:(kspan-R)
%for i=0:R:(kspan)
    for j=0:fsamp:(kspan-fsamp)
%    for j=0:fsamp:(kspan)
        counter = counter + 1;
        trajectory(counter,1) = j/FOV;
        trajectory(counter,2) = i/FOV;
    end
end

%Need to separate into odd and even cases, so that one point hits the
%origin.
nX = length(0:fsamp:(kspan-fsamp));
nY = length(0:R:(kspan-R));

kres = max(nX,nY);

xShift = mean(trajectory(:,1));
yShift = mean(trajectory(:,2));

trajectory(:,1) = trajectory(:,1) - xShift;
trajectory(:,2) = trajectory(:,2) - yShift;

if (mod(nX,1)==0) %if even
    trajectory(:,1) = trajectory(:,1) - fsamp/FOV/2;
    xShift = xShift + fsamp/FOV/2;
end

if (mod(nY,1)==0) %if even
    trajectory(:,2) = trajectory(:,2) - R/FOV/2;
    yShift = yShift + R/FOV/2;
end
