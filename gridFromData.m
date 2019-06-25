% after straightPolar

res = 256;
[trajectory, kres] = polarTraj(128,256,402,1,false);

tau = .0009;
ST = nuft_gg_init(2*pi*trajectory/res,[res res],8,2*res,tau); %Can add a tau variable

kpt = reshape(kdat,[numel(kdat) 1]);

x = nuft_gg_back(kpt,ST);

figure;imagesc(real(x))