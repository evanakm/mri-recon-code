
% Start with an analytic phantom
% Rasterize to the requested dimension
clear;

res = 256;
DefineBrain;
truth = RasterizePhantom(Brain,[res,res]);
%figure;imagesc(truth);colormap(pink)
%figure;imagesc(imresize(truth,[1024 1024]));colormap(pink)

DefineEdges;
edgeMap = RasterizePhantom(Edges,[res,res]);

sigma = .15;
Stage1=zeros(res,res);Stage1(edgeMap>0)=1;Stage1=imgaussfilt(Stage1,sigma);
EdgeReg=zeros(res,res);EdgeReg(Stage1 > 0)=1;
%figure;imagesc(EdgeReg)


DefineROI;
xyz = RasterizePhantom(Brain,[res,res]);




sigma = 10; alpha = .9;
%truth = RasterizePhantom(Brain,[res,res]);
%truth = imgaussfilt(phantom);

% Create k-space noise
noise = createCorrelatedNoise(sigma,alpha,res,res,res/2,res/2) + 1i*createCorrelatedNoise(sigma,alpha,res,res,res/2,res/2);
%noise = normrnd(0,sigma,res,res) + 1i*normrnd(0,sigma,res,res);
noNoise = zeros(size(noise));


% %There has to be a better way of doing this, but Matlab is being annoying.
%[indsX,indsY] = meshgrid(1:256,1:256);
% i1 = intersect(find(indsX>res*65/256),find(indsX<res*180/256));
% i1 = intersect(i1,find(indsY>res*60/256));
% i1 = intersect(i1,find(indsY<res*210/256));
% i1 = setdiff(1:res^2,i1);
% i2 = find(truth<0.001);
% mask = intersect(i1,i2);
% antimask = setdiff(1:res^2,mask);

%DefineROI;
%xyz = RasterizePhantom(Brain,[res res]);

mask =find(xyz == 0);
antimask = find(xyz > 0);
%antimask = 1:(256^2);

%SNR = @(image,truth) 10^(snr(abs(image(antimask)),abs(image(antimask))-abs(truth(antimask)))/20);
%SNR = @(image,truth) 10^(snr(abs(image(antimask)),abs(image(antimask))-abs(truth(antimask)))/20);

%test = ones([res res]);test(mask) = 0;figure;imagesc(test);

% Feed the rasterized and the noise through the model
% Baseline is fft --> addnoise --> ifft
X1 = truth; K1 = fftshift(fft2(truth)); Knoisy = K1+noise;
A1 = ifft2(K1); %Will be used as baseline for SNR calculations.
A2 = ifft2(Knoisy);
N1 = ifft2(noise);


figure;title('various plots');
subplot(2,2,1);imagesc(abs(X1));title('truth');
subplot(2,2,2);imagesc(abs(A1));title('Baseline Reconstruction');
subplot(2,2,3);imagesc(abs(A2));title('Noisy Reconstruction');
%subplot(2,2,3);imagesc(abs(X1)-abs(A1));title('noise');


% ----- OPTIONS FOR RECONSTRUCTION -------

counter = 0;

halfSpace = [];
halfSpace.run = true;

if(halfSpace.run)
    counter = counter + 1;
    
    direction = 'kx';
    Ahalf = halfKspace(Knoisy,direction);
    Nhalf = halfKspace(noise,direction);
    figure;subplot(2,2,1);imagesc(abs(X1));title('truth');colorbar;
    subplot(2,2,2);imagesc(abs(A2));title('reconstructed - full K-space');colorbar;
    subplot(2,2,3);imagesc(abs(Ahalf));title('reconstructed - half K-space');colorbar;
    subplot(2,2,4);imagesc(abs(X1)-abs(Ahalf));title('noise');

    sn = SNR(Ahalf,Ahalf-Nhalf,antimask);
    mse = roimse(Ahalf,A1,antimask);
    [ ssimAll, ssimROI ] = SSIM(Ahalf,X1,mask);
    ep = edgePenalty(Ahalf,X1,EdgeReg);
%    sn = SNR(Ahalf,X1);
%    mse = roimse(Ahalf,X1);
    
    eval(strcat('label',num2str(counter),' = ''Conjugate symmetry'';'));
    eval(strcat('snr',num2str(counter),' = sn;'));
    eval(strcat('mse',num2str(counter),' = mse;'));
    eval(strcat('ssimAll',num2str(counter),' = ssimAll;'));
    eval(strcat('ssimROI',num2str(counter),' = ssimROI;'));
    eval(strcat('ep',num2str(counter),' = ep;'));
    
end
    
homodyne = [];
homodyne.run = true;
homodyne.midband = 10;

if(homodyne.run)
    Astep = homodyneX(Knoisy,homodyne.midband,'step');
    Aramp = homodyneX(Knoisy,homodyne.midband,'ramp');
%    Astep = Astep';Aramp = Aramp';

    Nstep = homodyneX(noise,homodyne.midband,'step');
    Nramp = homodyneX(noise,homodyne.midband,'ramp');

    figure;subplot(2,3,1);imagesc(abs(X1));title('truth');colorbar;
    subplot(2,3,2);imagesc(abs(Astep));title('reconstructed - step');colorbar;
    subplot(2,3,3);imagesc(abs(Aramp));title('reconstructed - ramp');colorbar;
    
    subplot(2,3,5);imagesc(abs(X1)-abs(Astep));title('noise');colorbar;
    subplot(2,3,6);imagesc(abs(X1)-abs(Aramp));title('noise');colorbar;

    counter = counter + 1;
    sn = SNR(Astep,Astep-Nstep,antimask);
    mse = roimse(Astep,A1,antimask);
    [ ssimAll, ssimROI ] = SSIM(Astep,X1,mask);
    ep = edgePenalty(Astep,X1,EdgeReg);
    eval(strcat('label',num2str(counter),' = ''Homodyne with step'';'));
    eval(strcat('snr',num2str(counter),' = sn;'));
    eval(strcat('mse',num2str(counter),' = mse;'));
    eval(strcat('ssimAll',num2str(counter),' = ssimAll;'));
    eval(strcat('ssimROI',num2str(counter),' = ssimROI;'));
    eval(strcat('ep',num2str(counter),' = ep;'));
    
    counter = counter + 1;
    sn = SNR(Aramp,Aramp-Nramp,antimask);
    mse = roimse(Aramp,A1,antimask);
    [ ssimAll, ssimROI ] = SSIM(Aramp,X1,mask);
    ep = edgePenalty(Aramp,X1,EdgeReg);
    eval(strcat('label',num2str(counter),' = ''Homodyne with ramp'';'));
    eval(strcat('snr',num2str(counter),' = sn;'));
    eval(strcat('mse',num2str(counter),' = mse;'));
    eval(strcat('ssimAll',num2str(counter),' = ssimAll;'));
    eval(strcat('ssimROI',num2str(counter),' = ssimROI;'));
    eval(strcat('ep',num2str(counter),' = ep;'));
    
end

filteredBP = [];
filteredBP.runs = 4; % Can be a number between 0 and 4

filteredBP.run1.rStepSize = 1/2;
filteredBP.run1.angles = 402;
filteredBP.run1.filter = 'ram-lak';
filteredBP.run1.comparePlots = true;

filteredBP.run2.rStepSize = 1/2;
filteredBP.run2.angles = 201;
filteredBP.run2.filter = 'ram-lak';
filteredBP.run2.comparePlots = false;

filteredBP.run3.rStepSize = 1/2;
filteredBP.run3.angles = 100;
filteredBP.run3.filter = 'ram-lak';
filteredBP.run3.comparePlots = false;

filteredBP.run4.rStepSize = 1/2;
filteredBP.run4.angles = 50;
filteredBP.run4.filter = 'ram-lak';
filteredBP.run4.comparePlots = false;


    runs = filteredBP.runs;
    for i = 1:runs
        eval(strcat('Afbp',num2str(i),'=filteredBackProjection(X1,Knoisy,filteredBP.run',num2str(i),');'));
        eval(strcat('Afbp',num2str(i),'=imresize(Afbp',num2str(i),',[res res]);'));
        eval(strcat('[Afbp',num2str(i),',maxim]=normalize(Afbp',num2str(i),');'));

        eval(strcat('Nfbp',num2str(i),'=filteredBackProjection(X1,noise,filteredBP.run',num2str(i),');'));
        eval(strcat('Nfbp',num2str(i),'=imresize(Nfbp',num2str(i),',[res res]);'));
        eval(strcat('Nfbp',num2str(i),'=Nfbp',num2str(i),'/maxim;'));
    end

    figure;title('Filtered Back Projection');
    switch runs
        case 1
            subplot(2,2,1);imagesc(abs(X1));title('truth');colorbar;
            subplot(2,2,2);imagesc(abs(Afbp1));title('FBP - case 1');colorbar;
            subplot(2,2,4);imagesc(abs(X1)-abs(Afbp1));title('noise - case 1');colorbar;
            
            counter = counter + 1;
            sn = SNR(Afbp1,Afbp1-Nfbp1,antimask);
            mse = roimse(Afbp1,A1,antimask);
            [ ssimAll, ssimROI ] = SSIM(Afbp1,X1,mask);
            ep = edgePenalty(Afbp1,X1,EdgeReg);
            eval(strcat('label',num2str(counter),' = ''FBP scenario 1'';'));
            eval(strcat('snr',num2str(counter),' = sn;'));
            eval(strcat('mse',num2str(counter),' = mse;'));
            eval(strcat('ssimAll',num2str(counter),' = ssimAll;'));
            eval(strcat('ssimROI',num2str(counter),' = ssimROI;'));
            eval(strcat('ep',num2str(counter),' = ep;'));
    
        case 2
            subplot(2,3,1);imagesc(abs(X1));title('truth');colorbar;
            subplot(2,3,2);imagesc(abs(Afbp1));title('FBP - case 1');colorbar;
            subplot(2,3,5);imagesc(abs(X1)-abs(Afbp1));title('noise - case 1');colorbar;
            subplot(2,3,3);imagesc(abs(Afbp2));title('FBP - case 2');colorbar;
            subplot(2,3,6);imagesc(abs(X1)-abs(Afbp2));title('noise - case 2');colorbar;

            counter = counter + 1;
            sn = SNR(Afbp1,Afbp1-Nfbp1,antimask);
            mse = roimse(Afbp1,A1,antimask);
            [ ssimAll, ssimROI ] = SSIM(Afbp1,X1,mask);
            ep = edgePenalty(Afbp1,X1,EdgeReg);
            eval(strcat('label',num2str(counter),' = ''FBP scenario 1'';'));
            eval(strcat('snr',num2str(counter),' = sn;'));
            eval(strcat('mse',num2str(counter),' = mse;'));
            eval(strcat('ssimAll',num2str(counter),' = ssimAll;'));
            eval(strcat('ssimROI',num2str(counter),' = ssimROI;'));
            eval(strcat('ep',num2str(counter),' = ep;'));

            counter = counter + 1;
            sn = SNR(Afbp2,Afbp2-Nfbp2,antimask);
            mse = roimse(Afbp2,A1,antimask);
            [ ssimAll, ssimROI ] = SSIM(Afbp2,X1,mask);
            ep = edgePenalty(Afbp2,X1,EdgeReg);
            eval(strcat('label',num2str(counter),' = ''FBP scenario 2'';'));
            eval(strcat('snr',num2str(counter),' = sn;'));
            eval(strcat('mse',num2str(counter),' = mse;'));
            eval(strcat('ssimAll',num2str(counter),' = ssimAll;'));
            eval(strcat('ssimROI',num2str(counter),' = ssimROI;'));
            eval(strcat('ep',num2str(counter),' = ep;'));
    
        case 3
            subplot(2,4,1);imagesc(abs(X1));title('truth');colorbar;
            subplot(2,4,2);imagesc(abs(Afbp1));title('FBP - case 1');colorbar;
            subplot(2,4,6);imagesc(abs(X1)-abs(Afbp1));title('noise - case 1');colorbar;
            subplot(2,4,3);imagesc(abs(Afbp2));title('FBP - case 2');colorbar;
            subplot(2,4,7);imagesc(abs(X1)-abs(Afbp2));title('noise - case 2');colorbar;
            subplot(2,4,4);imagesc(abs(Afbp3));title('FBP - case 3');colorbar;
            subplot(2,4,8);imagesc(abs(X1)-abs(Afbp3));title('noise - case 3');colorbar;

            counter = counter + 1;
            sn = SNR(Afbp1,Afbp1-Nfbp1,antimask);
            mse = roimse(Afbp1,A1,antimask);
            [ ssimAll, ssimROI ] = SSIM(Afbp1,X1,mask);
            ep = edgePenalty(Afbp1,X1,EdgeReg);
            eval(strcat('label',num2str(counter),' = ''FBP scenario 1'';'));
            eval(strcat('snr',num2str(counter),' = sn;'));
            eval(strcat('mse',num2str(counter),' = mse;'));
            eval(strcat('ssimAll',num2str(counter),' = ssimAll;'));
            eval(strcat('ssimROI',num2str(counter),' = ssimROI;'));
            eval(strcat('ep',num2str(counter),' = ep;'));

            counter = counter + 1;
            sn = SNR(Afbp2,Afbp2-Nfbp2,antimask);
            mse = roimse(Afbp2,A1,antimask);
            [ ssimAll, ssimROI ] = SSIM(Afbp2,X1,mask);
            ep = edgePenalty(Afbp2,X1,EdgeReg);
            eval(strcat('label',num2str(counter),' = ''FBP scenario 2'';'));
            eval(strcat('snr',num2str(counter),' = sn;'));
            eval(strcat('mse',num2str(counter),' = mse;'));
            eval(strcat('ssimAll',num2str(counter),' = ssimAll;'));
            eval(strcat('ssimROI',num2str(counter),' = ssimROI;'));
            eval(strcat('ep',num2str(counter),' = ep;'));

            counter = counter + 1;
            sn = SNR(Afbp3,Afbp3-Nfbp3,antimask);
            mse = roimse(Afbp3,A1,antimask);
            [ ssimAll, ssimROI ] = SSIM(Afbp3,X1,mask);
            ep = edgePenalty(Afbp3,X1,EdgeReg);
            eval(strcat('label',num2str(counter),' = ''FBP scenario 3'';'));
            eval(strcat('snr',num2str(counter),' = sn;'));
            eval(strcat('mse',num2str(counter),' = mse;'));
            eval(strcat('ssimAll',num2str(counter),' = ssimAll;'));
            eval(strcat('ssimROI',num2str(counter),' = ssimROI;'));
            eval(strcat('ep',num2str(counter),' = ep;'));
            
        case 4
            subplot(4,3,1);imagesc(abs(X1));title('truth');colorbar;
            subplot(4,3,2);imagesc(abs(Afbp1));title('FBP - case 1');colorbar;
            subplot(4,3,5);imagesc(abs(X1)-abs(Afbp1));title('noise - case 1');colorbar;
            subplot(4,3,3);imagesc(abs(Afbp2));title('FBP - case 2');colorbar;
            subplot(4,3,6);imagesc(abs(X1)-abs(Afbp2));title('noise - case 2');colorbar;
            subplot(4,3,8);imagesc(abs(Afbp3));title('FBP - case 3');colorbar;
            subplot(4,3,11);imagesc(abs(X1)-abs(Afbp3));title('noise - case 3');colorbar;
            subplot(4,3,9);imagesc(abs(Afbp4));title('FBP - case 4');colorbar;
            subplot(4,3,12);imagesc(abs(X1)-abs(Afbp4));title('noise - case 4');colorbar;

            counter = counter + 1;
            sn = SNR(Afbp1,Afbp1-Nfbp1,antimask);
            mse = roimse(Afbp1,A1,antimask);
            [ ssimAll, ssimROI ] = SSIM(Afbp1,X1,mask);
            ep = edgePenalty(Afbp1,X1,EdgeReg);
            eval(strcat('label',num2str(counter),' = ''FBP scenario 1'';'));
            eval(strcat('snr',num2str(counter),' = sn;'));
            eval(strcat('mse',num2str(counter),' = mse;'));
            eval(strcat('ssimAll',num2str(counter),' = ssimAll;'));
            eval(strcat('ssimROI',num2str(counter),' = ssimROI;'));
            eval(strcat('ep',num2str(counter),' = ep;'));

            counter = counter + 1;
            sn = SNR(Afbp2,Afbp2-Nfbp2,antimask);
            mse = roimse(Afbp2,A1,antimask);
            [ ssimAll, ssimROI ] = SSIM(Afbp2,X1,mask);
            ep = edgePenalty(Afbp2,X1,EdgeReg);
            eval(strcat('label',num2str(counter),' = ''FBP scenario 2'';'));
            eval(strcat('snr',num2str(counter),' = sn;'));
            eval(strcat('mse',num2str(counter),' = mse;'));
            eval(strcat('ssimAll',num2str(counter),' = ssimAll;'));
            eval(strcat('ssimROI',num2str(counter),' = ssimROI;'));
            eval(strcat('ep',num2str(counter),' = ep;'));

            counter = counter + 1;
            sn = SNR(Afbp3,Afbp3-Nfbp3,antimask);
            mse = roimse(Afbp3,A1,antimask);
            [ ssimAll, ssimROI ] = SSIM(Afbp3,X1,mask);
            ep = edgePenalty(Afbp3,X1,EdgeReg);
            eval(strcat('label',num2str(counter),' = ''FBP scenario 3'';'));
            eval(strcat('snr',num2str(counter),' = sn;'));
            eval(strcat('mse',num2str(counter),' = mse;'));
            eval(strcat('ssimAll',num2str(counter),' = ssimAll;'));
            eval(strcat('ssimROI',num2str(counter),' = ssimROI;'));
            eval(strcat('ep',num2str(counter),' = ep;'));
            
            counter = counter + 1;
            sn = SNR(Afbp4,Afbp4-Nfbp4,antimask);
            mse = roimse(Afbp4,A1,antimask);
            [ ssimAll, ssimROI ] = SSIM(Afbp4,X1,mask);
            ep = edgePenalty(Afbp4,X1,EdgeReg);
            eval(strcat('label',num2str(counter),' = ''FBP scenario 4'';'));
            eval(strcat('snr',num2str(counter),' = sn;'));
            eval(strcat('mse',num2str(counter),' = mse;'));
            eval(strcat('ssimAll',num2str(counter),' = ssimAll;'));
            eval(strcat('ssimROI',num2str(counter),' = ssimROI;'));
            eval(strcat('ep',num2str(counter),' = ep;'));
            
        otherwise
                 
    end


gridding = [];
gridding.runs = 2; % between 0 and 4;

% some notes: For a cartesian trajectory,
%                   par1 = step size in kX direction, internally scaled so
%                          that 1 is the Nyquist
%                   par2 = step size in kY direction
%             For a polar trajectory
%                   par1 = number of steps in the R direction
%                   par2 = number of steps in the azimuthal direction
%             For a spiral trajectory
%                   par1 = number of steps per leaf
%                   par2 = number of revolutions per leaf
%                   par3 = number of interleaved paths
%             For a custom trajectory
%                   par1 = an array of [kX,kY] pairs

%gridding.run1.interp = 'linear';
%gridding.run1.traj = 'cartesian';
%gridding.run1.tau = .15; %Experimentally, .15 seems to give good images
%gridding.run1.res = res;
%gridding.run1.kReach = res/2;
%gridding.run1.nearest = 8; %number of nearest neighbors in the internal gaussian gridding algorithm. Tradeoff between speed and accuracy.
%gridding.run1.over = 3;
%gridding.run1.par1 = 1;
%gridding.run1.par2 = 1;

gridding.run1.interp = 'linear';
gridding.run1.traj = 'polar';
gridding.run1.tau = .15; %Experimentally, .15 seems to give good images
gridding.run1.res = res;
gridding.run1.kReach = res/2;
gridding.run1.nearest = 3; %number of nearest neighbors in the internal gaussian gridding algorithm. Tradeoff between speed and accuracy.
gridding.run1.over = 3;
gridding.run1.par1 = .8*2*1.3*res;
%gridding.run2.par1 = 2*res;
gridding.run1.par2 = 402;


gridding.run2.interp = 'linear';
gridding.run2.traj = 'polar';
gridding.run2.tau = .15; %Experimentally, .15 seems to give good images
gridding.run2.res = res;
gridding.run2.kReach = res/2;
gridding.run2.nearest = 8; %number of nearest neighbors in the internal gaussian gridding algorithm. Tradeoff between speed and accuracy.
gridding.run2.over = 3;
gridding.run2.par1 = 2*1.3*res;
%gridding.run2.par1 = 2*res;
gridding.run2.par2 = 402;

gridding.run3.interp = 'gaussian';
gridding.run3.traj = 'cartesian';
gridding.run3.tau = .15; %Experimentally, .15 seems to give good images
gridding.run3.res = res;
gridding.run3.kReach = res/2;
gridding.run3.nearest = 8; %number of nearest neighbors in the internal gaussian gridding algorithm. Tradeoff between speed and accuracy.
gridding.run3.over = 3;
gridding.run3.par1 = 1;
gridding.run3.par2 = 1;

gridding.run4.interp = 'gaussian';
gridding.run4.traj = 'polar';
gridding.run4.tau = .15; %Experimentally, .15 seems to give good images
gridding.run4.res = res;
gridding.run4.kReach = res/2;
gridding.run4.over = 3;
gridding.run4.nearest = 8; %number of nearest neighbors in the internal gaussian gridding algorithm. Tradeoff between speed and accuracy.
gridding.run4.par1 = 2*1.3*res;
gridding.run4.par2 = 2*402;

% gridding.run1.interp = 'linear';
% gridding.run1.traj = 'spiral';
% gridding.run1.tau = .15; %Experimentally, .15 seems to give good images
% gridding.run1.res = res;
% gridding.run1.kReach = res/2;
% gridding.run1.nearest = 8; %number of nearest neighbors in the internal gaussian gridding algorithm. Tradeoff between speed and accuracy.
% gridding.run1.over = 3;
% gridding.run1.par1 = 800;
% gridding.run1.par2 = 100;
% gridding.run1.par3 = 10000;
% 
% gridding.run2.interp = 'linear';
% gridding.run2.traj = 'spiral';
% gridding.run2.tau = .15; %Experimentally, .15 seems to give good images
% gridding.run2.res = res;
% gridding.run2.kReach = res/2;
% gridding.run2.nearest = 8; %number of nearest neighbors in the internal gaussian gridding algorithm. Tradeoff between speed and accuracy.
% gridding.run2.over = 3;
% gridding.run2.par1 = 80000;
% gridding.run2.par2 = 10000;
% gridding.run2.par3 = 1;
% 
% gridding.run3.interp = 'gaussian';
% gridding.run3.traj = 'spiral';
% gridding.run3.tau = .15; %Experimentally, .15 seems to give good images
% gridding.run3.res = res;
% gridding.run3.kReach = res/2;
% gridding.run3.nearest = 8; %number of nearest neighbors in the internal gaussian gridding algorithm. Tradeoff between speed and accuracy.
% gridding.run3.over = 3;
% gridding.run3.par1 = 800;
% gridding.run3.par2 = 100;
% gridding.run3.par3 = 100;
% 
% gridding.run4.interp = 'gaussian';
% gridding.run4.traj = 'spiral';
% gridding.run4.tau = .15; %Experimentally, .15 seems to give good images
% gridding.run4.res = res;
% gridding.run4.kReach = res/2;
% gridding.run4.nearest = 8; %number of nearest neighbors in the internal gaussian gridding algorithm. Tradeoff between speed and accuracy.
% gridding.run4.over = 3;
% gridding.run4.par1 = 80000;
% gridding.run4.par2 = 10000;
% gridding.run4.par3 = 1;

if(gridding.runs>0)
    runs = gridding.runs;
    for i = 1:runs
%        eval(strcat('Agrid',num2str(i),'=griddingAlgo(Knoisy,gridding.run',num2str(i),');'));
        eval(strcat('Agrid',num2str(i),'=griddingAlgo(K1,gridding.run',num2str(i),');'));
        disp(strcat('Gridding run number',{' '},num2str(i),' complete.'));
        eval(strcat('[Agrid',num2str(i),',maxim]=normalize(Agrid',num2str(i),');'));
        eval(strcat('Ngrid',num2str(i),'=griddingAlgo(noise,gridding.run',num2str(i),')/maxim;'));

    end

    figure;title('Gridding Algorithm');
    switch runs
        case 1
            subplot(2,2,1);imagesc(abs(X1));title('truth');colorbar;
            subplot(2,2,2);imagesc(abs(Agrid1));title('Gridding - case 1');colorbar;
            subplot(2,2,4);imagesc(abs(X1)-abs(Agrid1));title('noise - case 1');colorbar;

            counter = counter + 1;
            sn = SNR(Agrid1,Agrid1-Ngrid1,antimask);
            mse = roimse(Agrid1,A1,antimask);
            [ ssimAll, ssimROI ] = SSIM(Agrid1,X1,mask);
            ep = edgePenalty(Agrid1,X1,EdgeReg);
            eval(strcat('label',num2str(counter),' = ''Gridding scenario 1'';'));
            eval(strcat('snr',num2str(counter),' = sn;'));
            eval(strcat('mse',num2str(counter),' = mse;'));
            eval(strcat('ssimAll',num2str(counter),' = ssimAll;'));
            eval(strcat('ssimROI',num2str(counter),' = ssimROI;'));
            eval(strcat('ep',num2str(counter),' = ep;'));
            
        case 2
            subplot(2,3,1);imagesc(abs(X1));title('truth');colorbar;
            subplot(2,3,2);imagesc(abs(Agrid1));title('Gridding - case 1');colorbar;
            subplot(2,3,5);imagesc(abs(X1)-abs(Agrid1));title('noise - case 1');colorbar;
            subplot(2,3,3);imagesc(abs(Agrid2));title('Gridding - case 2');colorbar;
            subplot(2,3,6);imagesc(abs(X1)-abs(Agrid2));title('noise - case 2');colorbar;

            counter = counter + 1;
            sn = SNR(Agrid1,Agrid1-Ngrid1,antimask);
            mse = roimse(Agrid1,A1,antimask);
            [ ssimAll, ssimROI ] = SSIM(Agrid1,X1,mask);
            ep = edgePenalty(Agrid1,X1,EdgeReg);
            eval(strcat('label',num2str(counter),' = ''Gridding scenario 1'';'));
            eval(strcat('snr',num2str(counter),' = sn;'));
            eval(strcat('mse',num2str(counter),' = mse;'));
            eval(strcat('ssimAll',num2str(counter),' = ssimAll;'));
            eval(strcat('ssimROI',num2str(counter),' = ssimROI;'));
            eval(strcat('ep',num2str(counter),' = ep;'));

            counter = counter + 1;
            sn = SNR(Agrid2,Agrid2-Ngrid2,antimask);
            mse = roimse(Agrid2,A1,antimask);
            [ ssimAll, ssimROI ] = SSIM(Agrid2,X1,mask);
            ep = edgePenalty(Agrid2,X1,EdgeReg);
            eval(strcat('label',num2str(counter),' = ''Gridding scenario 2'';'));
            eval(strcat('snr',num2str(counter),' = sn;'));
            eval(strcat('mse',num2str(counter),' = mse;'));
            eval(strcat('ssimAll',num2str(counter),' = ssimAll;'));
            eval(strcat('ssimROI',num2str(counter),' = ssimROI;'));
            eval(strcat('ep',num2str(counter),' = ep;'));
            
        case 3
            subplot(2,4,1);imagesc(abs(X1));title('truth');colorbar;
            subplot(2,4,2);imagesc(abs(Agrid1));title('Gridding - case 1');colorbar;
            subplot(2,4,6);imagesc(abs(X1)-abs(Agrid1));title('noise - case 1');colorbar;
            subplot(2,4,3);imagesc(abs(Agrid2));title('Gridding - case 2');colorbar;
            subplot(2,4,7);imagesc(abs(X1)-abs(Agrid2));title('noise - case 2');colorbar;
            subplot(2,4,4);imagesc(abs(Agrid3));title('Gridding - case 3');colorbar;
            subplot(2,4,8);imagesc(abs(X1)-abs(Agrid3));title('noise - case 3');colorbar;
            
            counter = counter + 1;
            sn = SNR(Agrid1,Agrid1-Ngrid1,antimask);
            mse = roimse(Agrid1,A1,antimask);
            [ ssimAll, ssimROI ] = SSIM(Agrid1,X1,mask);
            ep = edgePenalty(Agrid1,X1,EdgeReg);
            eval(strcat('label',num2str(counter),' = ''Gridding scenario 1'';'));
            eval(strcat('snr',num2str(counter),' = sn;'));
            eval(strcat('mse',num2str(counter),' = mse;'));
            eval(strcat('ssimAll',num2str(counter),' = ssimAll;'));
            eval(strcat('ssimROI',num2str(counter),' = ssimROI;'));
            eval(strcat('ep',num2str(counter),' = ep;'));

            counter = counter + 1;
            sn = SNR(Agrid2,Agrid2-Ngrid2,antimask);
            mse = roimse(Agrid2,A1,antimask);
            [ ssimAll, ssimROI ] = SSIM(Agrid2,X1,mask);
            ep = edgePenalty(Agrid2,X1,EdgeReg);
            eval(strcat('label',num2str(counter),' = ''Gridding scenario 2'';'));
            eval(strcat('snr',num2str(counter),' = sn;'));
            eval(strcat('mse',num2str(counter),' = mse;'));
            eval(strcat('ssimAll',num2str(counter),' = ssimAll;'));
            eval(strcat('ssimROI',num2str(counter),' = ssimROI;'));
            eval(strcat('ep',num2str(counter),' = ep;'));

            counter = counter + 1;
            sn = SNR(Agrid3,Agrid3-Ngrid3,antimask);
            mse = roimse(Agrid3,A1,antimask);
            [ ssimAll, ssimROI ] = SSIM(Agrid3,X1,mask);
            ep = edgePenalty(Agrid3,X1,EdgeReg);
            eval(strcat('label',num2str(counter),' = ''Gridding scenario 3'';'));
            eval(strcat('snr',num2str(counter),' = sn;'));
            eval(strcat('mse',num2str(counter),' = mse;'));
            eval(strcat('ssimAll',num2str(counter),' = ssimAll;'));
            eval(strcat('ssimROI',num2str(counter),' = ssimROI;'));
            eval(strcat('ep',num2str(counter),' = ep;'));
            
        case 4
            subplot(4,3,1);imagesc(abs(X1));title('truth');colorbar;
            subplot(4,3,2);imagesc(abs(Agrid1));title('Gridding - case 1');colorbar;
            subplot(4,3,5);imagesc(abs(X1)-abs(Agrid1));title('noise - case 1');colorbar;
            subplot(4,3,3);imagesc(abs(Agrid2));title('Gridding - case 2');colorbar;
            subplot(4,3,6);imagesc(abs(X1)-abs(Agrid2));title('noise - case 2');colorbar;
            subplot(4,3,8);imagesc(abs(Agrid3));title('Gridding - case 3');colorbar;
            subplot(4,3,11);imagesc(abs(X1)-abs(Agrid3));title('noise - case 3');colorbar;
            subplot(4,3,9);imagesc(abs(Agrid4));title('Gridding - case 4');colorbar;
            subplot(4,3,12);imagesc(abs(X1)-abs(Agrid4));title('noise - case 4');colorbar;

            counter = counter + 1;
            sn = SNR(Agrid1,Agrid1-Ngrid1,antimask);
            mse = roimse(Agrid1,A1,antimask);
            [ ssimAll, ssimROI ] = SSIM(Agrid1,X1,mask);
            ep = edgePenalty(Agrid1,X1,EdgeReg);
            eval(strcat('label',num2str(counter),' = ''Gridding scenario 1'';'));
            eval(strcat('snr',num2str(counter),' = sn;'));
            eval(strcat('mse',num2str(counter),' = mse;'));
            eval(strcat('ssimAll',num2str(counter),' = ssimAll;'));
            eval(strcat('ssimROI',num2str(counter),' = ssimROI;'));
            eval(strcat('ep',num2str(counter),' = ep;'));

            counter = counter + 1;
            sn = SNR(Agrid2,Agrid2-Ngrid2,antimask);
            mse = roimse(Agrid2,A1,antimask);
            [ ssimAll, ssimROI ] = SSIM(Agrid2,X1,mask);
            ep = edgePenalty(Agrid2,X1,EdgeReg);
            eval(strcat('label',num2str(counter),' = ''Gridding scenario 2'';'));
            eval(strcat('snr',num2str(counter),' = sn;'));
            eval(strcat('mse',num2str(counter),' = mse;'));
            eval(strcat('ssimAll',num2str(counter),' = ssimAll;'));
            eval(strcat('ssimROI',num2str(counter),' = ssimROI;'));
            eval(strcat('ep',num2str(counter),' = ep;'));

            counter = counter + 1;
            sn = SNR(Agrid3,Agrid3-Ngrid3,antimask);
            mse = roimse(Agrid3,A1,antimask);
            [ ssimAll, ssimROI ] = SSIM(Agrid3,X1,mask);
            ep = edgePenalty(Agrid3,X1,EdgeReg);
            eval(strcat('label',num2str(counter),' = ''Gridding scenario 3'';'));
            eval(strcat('snr',num2str(counter),' = sn;'));
            eval(strcat('mse',num2str(counter),' = mse;'));
            eval(strcat('ssimAll',num2str(counter),' = ssimAll;'));
            eval(strcat('ssimROI',num2str(counter),' = ssimROI;'));
            eval(strcat('ep',num2str(counter),' = ep;'));
            
            counter = counter + 1;
            sn = SNR(Agrid4,Agrid4-Ngrid4,antimask);
            mse = roimse(Agrid4,A1,antimask);
            [ ssimAll, ssimROI ] = SSIM(Agrid4,X1,mask);
            ep = edgePenalty(Agrid4,X1,EdgeReg);
            eval(strcat('label',num2str(counter),' = ''Gridding scenario 4'';'));
            eval(strcat('snr',num2str(counter),' = sn;'));
            eval(strcat('mse',num2str(counter),' = mse;'));
            eval(strcat('ssimAll',num2str(counter),' = ssimAll;'));
            eval(strcat('ssimROI',num2str(counter),' = ssimROI;'));
            eval(strcat('ep',num2str(counter),' = ep;'));
                        
    
    end
end

snr = cell(counter+2,7);
snr{1,1} = 'RUN';
snr{1,2} = 'PSNR';
snr{1,3} = 'MSE';
snr{1,4} = 'SSIM - all';
snr{1,5} = 'SSIM - ROI';
snr{1,6} = 'Edge Penalty';

snr{2,1} = 'Baseline';
base = SNR(A2,A1,antimask);
snr{2,2} = base.snr;
snr{2,3} = roimse(A2,A1,antimask);
[all, roiOnly] = SSIM(A2,X1,mask);
snr{2,4} = all;
snr{2,5} = roiOnly;
snr{2,6} = edgePenalty(A2,X1,EdgeReg);

for i = 1:counter
    eval(strcat('snr{',num2str(i+2),',1} = label',num2str(i),';'));
    eval(strcat('snr{',num2str(i+2),',2} = snr',num2str(i),'.snr;'));
    eval(strcat('snr{',num2str(i+2),',3} = mse',num2str(i),';'));
    eval(strcat('snr{',num2str(i+2),',4} = ssimAll',num2str(i),';'));
    eval(strcat('snr{',num2str(i+2),',5} = ssimROI',num2str(i),';'));
    eval(strcat('snr{',num2str(i+2),',6} = ep',num2str(i),';'));
    eval(strcat('snr{',num2str(i+2),',7} = ssimROI',num2str(i),'*(1-.5*(1-ep',num2str(i),'));'));
%    eval(strcat('snr{',num2str(i+2),',6} = snr',num2str(i),'.signal;'));
%    eval(strcat('snr{',num2str(i+2),',7} = snr',num2str(i),'.noise;'));
end

snr
% PRINT OUT SNR
% PRINT OUT MSE
% PRINT OUT SSIM


%Astep(mask)=0;figure;imagesc(abs(abs(truth)-abs(Astep)));colormap(summer);title('Homodyne (step)')

