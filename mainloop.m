clear;

kernel = @(x,a) max(1 - abs(x/a),0);

res = 256;

%h = fspecial('disk',res/16);

truth = imresize(imgaussfilt(phantom),[res res]);
%truth = fftshift(truth);
base = ones([res res]);

kspace = fftshift(fftn(truth));
kbase = fftshift(fftn(base));
%kspace = fftn(truth);
%kbase = fftn(base);

[indsX,indsY] = meshgrid(1:256,1:256);

%There has to be a better way of doing this, but Matlab is being annoying.

i1 = intersect(find(indsX>65),find(indsX<180));
i1 = intersect(i1,find(indsY>60));
i1 = intersect(i1,find(indsY<210));
i1 = setdiff(1:256^2,i1);
i2 = find(truth<0.001);
mask = intersect(i1,i2);
test = ones([res res]);test(mask) = 0;%figure;imagesc(test);

%Now let's zero-pad
xpad = 0;
ypad = 0;

rows = size(kspace,1);
cols = size(kspace,2);

kspace = [zeros(rows,xpad), kspace, zeros(rows,xpad)];
kspace = [zeros(ypad,cols + 2*xpad); kspace; zeros(ypad,cols + 2*xpad)];

kbase = [zeros(rows,xpad), kbase, zeros(rows,xpad)];
kbase = [zeros(ypad,cols + 2*xpad); kbase; zeros(ypad,cols + 2*xpad)];

%res = FOV/256; %not sure if needed


%---------------------------------------------------

FOV = 1; %.28 in the code, but let's normalize.
kReach = res/2;
%rSteps = 400;
%kReach = 2*res;

rSteps = 1.5*res;thSteps = 400;
[trajectory, kres] = polarTraj(kReach,rSteps,thSteps,FOV);
[kPTS, kernK, kerntraj, gridsize] = triInterp(trajectory,res,FOV,kspace);
ST = nuft_gg_init(2*pi*trajectory/res,[res res],8,5*kres); %Can add a tau variable
x = nuft_gg_back(kPTS,ST);
%b = nuft_gg_back(basePTS,ST);
ST2 = nuft_gg_init(kerntraj,[res res],8,5*gridsize);
k = nuft_gg_back(kernK,ST2);

% ------------------------------- ALL TESTS FROM HERE ON
% -----------------------

figure;subplot(2,2,1);imagesc(abs(x)./abs(k));
%subplot(2,2,2);imagesc(abs(k));

fsamp = .25; R = .25;
[trajectory, kres] = cartesianTraj(kReach,fsamp,R,FOV);
[kPTS, kernK, kerntraj, gridsize] = triInterp(trajectory,res,FOV,kspace);
ST = nuft_gg_init(2*pi*trajectory/res,[res res],8,5*kres); %Can add a tau variable
x = nuft_gg_back(kPTS,ST);
ST2 = nuft_gg_init(kerntraj,[res res],8,5*gridsize);
k = nuft_gg_back(kernK,ST2);
subplot(2,2,2);imagesc(abs(x)./abs(k));

numSteps = 800;numLoops = 100;numLeaves = 100;
[trajectory, kres] = spiralTraj(kReach,numSteps,numLoops,numLeaves);
[kPTS, kernK, kerntraj, gridsize] = triInterp(trajectory,res,FOV,kspace);
ST = nuft_gg_init(2*pi*trajectory/res,[res res],8,3*res); %Can add a tau variable
x = nuft_gg_back(kPTS,ST);
ST2 = nuft_gg_init(kerntraj,[res res],8,5*gridsize);
k = nuft_gg_back(kernK,ST2);
subplot(2,2,3);imagesc(abs(x)./abs(k));

subplot(2,2,4);imagesc(abs(k));


bottom = -res/2; top = (res/2)-1;
[kXgrid,kYgrid] = meshgrid((bottom:top)/FOV,(bottom:top)/FOV);

tau = .15;

%rSteps = 1.5*res;thSteps = 400;
rSteps = 2*res;thSteps = 100;
[trajectory, kres] = polarTraj(kReach,rSteps,thSteps,FOV);
[kPTS, kernK, kerntraj, gridsize] = gaussianInterp(trajectory,kspace,kXgrid,kYgrid,tau,res);
ST = nuft_gg_init(2*pi*trajectory/res,[res res],8,5*kres); %Can add a tau variable
x = nuft_gg_back(kPTS,ST);
%b = nuft_gg_back(basePTS,ST);
ST2 = nuft_gg_init(kerntraj,[res res],8,5*gridsize);
k = nuft_gg_back(kernK,ST2);

figure; subplot(2,2,1); imagesc(abs(x)./abs(k));

fsamp = 1; R = 1;
[trajectory, kres] = cartesianTraj(kReach,fsamp,R,FOV);
[kPTS, kernK, kerntraj, gridsize] = gaussianInterp(trajectory,kspace,kXgrid,kYgrid,tau,res);
ST = nuft_gg_init(2*pi*trajectory/res,[res res],8,5*kres); %Can add a tau variable
x = nuft_gg_back(kPTS,ST);
ST2 = nuft_gg_init(kerntraj,[res res],8,5*gridsize);
k = nuft_gg_back(kernK,ST2);
subplot(2,2,2);imagesc(abs(x)./abs(k));


numSteps = 800;numLoops = 100;numLeaves = 100;
[trajectory, kres] = spiralTraj(kReach,numSteps,numLoops,numLeaves);
[kPTS, kernK, kerntraj, gridsize] = gaussianInterp(trajectory,kspace,kXgrid,kYgrid,tau,res);
ST = nuft_gg_init(2*pi*trajectory/res,[res res],8,3*res); %Can add a tau variable
x = nuft_gg_back(kPTS,ST);
ST2 = nuft_gg_init(kerntraj,[res res],8,5*gridsize);
k = nuft_gg_back(kernK,ST2);
subplot(2,2,3);imagesc(abs(x)./abs(k));


numSteps = 120000;numLoops = 250;numLeaves = 1;
[trajectory, kres] = spiralTraj(kReach,numSteps,numLoops,numLeaves);
[kPTS, kernK, kerntraj, gridsize] = gaussianInterp(trajectory,kspace,kXgrid,kYgrid,tau,res);
ST = nuft_gg_init(2*pi*trajectory/res,[res res],8,3*res); %Can add a tau variable
x = nuft_gg_back(kPTS,ST);
ST2 = nuft_gg_init(kerntraj,[res res],8,5*gridsize);
k = nuft_gg_back(kernK,ST2);
subplot(2,2,4);imagesc(abs(x)./abs(k));


%m = abs(x)';
%m(mask)=0;
%subplot(2,2,4);imagesc(abs(m)./abs(k));


%subplot(2,2,2); imagesc(abs(x));

%subplot(2,2,4); imagesc(abs(k));


% x = nuft_gg_back(kPTS,ST);
% %b = nuft_gg_back(basePTS,ST);
% ST2 = nuft_gg_init(kerntraj,[res res],8,5*gridsize);
% k = nuft_gg_back(kernK,ST2);
% 
% %There has to be a better way of doing this, but Matlab is being annoying.
% 
% i1 = intersect(find(indsX>65),find(indsX<180));
% i1 = intersect(i1,find(indsY>60));
% i1 = intersect(i1,find(indsY<210));
% i1 = setdiff(1:256^2,i1);
% i2 = find(truth<0.001);
% mask = intersect(i1,i2);
% test = ones([res res]);test(mask) = 0;%figure;imagesc(test);
% x(mask) = 0;
% 
% subplot(2,2,3);imagesc(abs(x)./abs(k));
% subplot(2,2,4);imagesc(abs(k));


%     %x=fftshift(x);
%     %b=fftshift(b);
% 
%     XX = linspace(0,1,res);
%     YY = linspace(0,1,res);
%     [meshX, meshY] = meshgrid(XX,YY);
% 
%     %XX = linspace(0,1,res);
%     %YY = linspace(0,1,res);
% 
%     %Need to cap it
%     %pic = abs(x)./abs(s2)./abs(s1);
%     %pic(mask) = 0;
% 
%     %maxpic = max(pic(:));
% 
%     %pic = min(pic,.9*maxpic);
% 
%     m = x;
%     m(mask) = 0;
% 
%     q = abs(x)./abs(b);
%     M = max(q(:));
%     q = min(q,.2*M);
% 
%     %figure;imagesc(pic);
%     figure;
%     subplot(2,2,1);imagesc(abs(x));
%     subplot(2,2,2);imagesc(abs(b));
%     subplot(2,2,3);imagesc(abs(x)./abs(b));
%     subplot(2,2,4);imagesc(abs(q));
% 
%     figure;
%     subplot(2,2,1);imagesc(abs(x));
%     subplot(2,2,2);imagesc(abs(b));
%     subplot(2,2,3);imagesc(abs(m));
%     subplot(2,2,4);imagesc(abs(m)./abs(b));
% 
%     factor = 4;
%     XX = linspace(-0.5,0.5,res);
%     YY = linspace(-0.5,0.5,res);
%     [meshX, meshY] =  meshgrid(XX,YY);
%     meshX = meshX/factor;
%     meshY = meshY/factor;
%     s1 = sinc(meshX/FOV).^2;
%     s2 = sinc(meshY/FOV).^2;
% 
%     s = s1 .* s2;
%     s2 = abs(x)./abs(s);
%     s2(mask) = 0;
% 
%     q2 = abs(x)./abs(s);
%     q2(mask) = 0;
% 
%     figure;
%     subplot(2,2,1);imagesc(abs(x));
%     subplot(2,2,2);imagesc(abs(s));colorbar;
%     subplot(2,2,3);imagesc(abs(x)./abs(s));
%     %subplot(2,2,4);imagesc(abs(q2));
% 
%     % p = abs(x)./abs(s1.*s2);
%     % M = max(p(:));
%     % p = min(p,.0005*M);
%     %
%     % figure;
%     % subplot(2,2,1);imagesc(abs(x));
%     % subplot(2,2,2);imagesc(abs(s1.*s2));
%     % subplot(2,2,3);imagesc(abs(x)./abs(s1.*s2));
%     % subplot(2,2,4);imagesc(abs(p));
%     %
%     %
%     % o = abs(s1.*s2);
%     % M = max(o(:));
%     % o = max(o,.005*M);
%     %
%     % o2 = abs(x)./abs(o);
%     % M2 = max(o2(:));
%     % o2 = min(o2,.005*M2);
%     %
%     % figure;
%     % subplot(2,2,1);imagesc(abs(x));
%     % subplot(2,2,2);imagesc(abs(o));
%     % subplot(2,2,3);imagesc(abs(x)./abs(o));
%     % subplot(2,2,4);imagesc(abs(o2));
%     %
% 
%     %subplot(2,2,3);imagesc(abs(q));
% 
%     %subplot(2,2,4);imagesc(abs(m));
% 
% 
% 
% 

