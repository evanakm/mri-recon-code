function [image] = griddingTriangular(kspace,trajectory,res,kres,nearest,FOV)

% Do not need to pad k-space in the current iteration of the code.

[kPTS, kernK, kerntraj, gridsize] = triInterp(trajectory,res,FOV,kspace);

ST = nuft_gg_init(2*pi*trajectory/res,[res res],nearest,3*res); %Can add a tau variable
x = nuft_gg_back(kPTS,ST);

% Need to invert the convolutional kernel as well.
ST2 = nuft_gg_init(kerntraj,[res res],nearest,5*gridsize);
k = nuft_gg_back(kernK,ST2);

image = abs(x)./abs(k);