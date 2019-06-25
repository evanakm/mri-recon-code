clear
fid = fopen('spgr3d_256x_256y_160z_1e_12c.raw','r');
hdr=157276;
Nx=256; Ny=256; Nz=160; Nc=12;

bpp=2; bpL=Nx*2*bpp;
bpS=bpL*(Ny+1);bpC=bpS*Nz;

slices = zeros(Nz,Nx,Ny);
sliceNoise = zeros(Nz,Nx,Ny);
sliceFilter = zeros(Nz,Nx,Ny);
sliceBoth = zeros(Nz,Nx,Ny);



%filter = zeros(Nx,Ny);
%[X,Y] = meshgrid(exp(1i*(-128:127)*pi/256),exp(1i*(-128:127)*pi/256));
%Z = X.*Y;
%[X,Y] = meshgrid(hamming(256),hamming(256));
[X,Y] = meshgrid(-128:127,-128:127);
Z = sqrt(X.^2+Y.^2);
Z = 0.5*max(Z(:)) - Z;
filt = max(Z,0);

%filt = Z/15;

step = 4;
inds = 1:step:255;

sliceUnder = zeros(Nz,Nx/step,Ny/step);
sliceUnderNoise = zeros(Nz,Nx/step,Ny/step);

%Saving k-space for debugging purposes
kspace = zeros(Nz,Nx,Ny);

for z=1:Nz
    for c=1:Nc
        fseek(fid,hdr+(c-1)*bpC+bpS*(z-1)+1*bpL,'bof');
        thisSlice = fread(fid,[2*Nx,Ny],'short');
        
        kdat = thisSlice(1:2:end-1,:) + 1i*thisSlice(2:2:end,:);
        
        kspace(z,:,:) = kspace(z,:,:) + reshape(kdat,size(kspace(z,:,:)));
        
        slices(z,:,:) = slices(z,:,:) + reshape(ifft2(kdat),[1 Nx Ny]);

        %Add gaussian noise
        kdat1 = kdat + normrnd(0,1,[Nx,Ny])+1i*normrnd(0,1,[Nx,Ny]);
        kdat2 = kdat .* filt;
        kdat3 = kdat1 .* filt;
%        kdat2 = kdat(inds,inds);
%        kdat3 = kdat1(inds,inds);



%        kdat = kdat .* exp(1i*normrnd(pi,2*pi,[Nx,Ny]));
%        runix = reshape(-128 + 255*rand(256*256,1),[256 256]);
%        kdat = kdat .* exp(1i*pi*(X+Y)*(1/256).*runix);
%        kdat = kdat .* filt;
        sliceNoise(z,:,:) = sliceNoise(z,:,:) + reshape(ifft2(kdat1),[1 Nx Ny]);
        sliceFilter(z,:,:) = sliceFilter(z,:,:) + reshape(ifft2(kdat2),[1 Nx Ny]);
        sliceBoth(z,:,:) = sliceBoth(z,:,:) + reshape(ifft2(kdat3),[1 Nx Ny]);
%        sliceUnder(z,:,:) = sliceUnder(z,:,:) + reshape(ifft2(kdat2),[1,Nx/step,Ny/step]);
%        sliceUnderNoise(z,:,:) = sliceUnderNoise(z,:,:) + reshape(ifft2(kdat3),[1,Nx/step,Ny/step]);


    end
end

fclose(fid);

figure
slice = 40;

F1 = fftshift(reshape(abs(slices(slice,:,:)),[Nx,Ny]));
F2 = fftshift(reshape(abs(sliceNoise(slice,:,:)),[Nx,Ny]));
F3 = fftshift(reshape(abs(sliceFilter(slice,:,:)),[Nx,Ny]));
F4 = fftshift(reshape(abs(sliceBoth(slice,:,:)),[Nx,Ny]));

subplot(2,2,1);
%imshow(fftshift(reshape(abs(slices(slice,:,:)),[Nx,Ny])));colorbar;
imagesc(F1);colorbar;
title('Original');
subplot(2,2,2);
imagesc(F2);colorbar;
title('With Noise');
subplot(2,2,3);
imagesc(F3);colorbar;
title('With Filter');
subplot(2,2,4);
imagesc(F4);colorbar;
title('With Both');
%subplot(2,2,3);
%imagesc(fftshift(reshape(abs(sliceUnder(slice,:,:)),[Nx/step,Ny/step])));colorbar;
%title('Undersampled');
%subplot(2,2,4);
%imagesc(fftshift(reshape(abs(sliceUnderNoise(slice,:,:)),[Nx/step,Ny/step])));colorbar;
%title('Undersampled With Noise');
%--------------

thr = .025;

SNR = 10^(snr(F1,F2-F1)/20);
roimat = F1(F1>thr); roinoise = F2(F1>thr) - roimat;

ROIsnr = 10^(snr(roimat,roinoise)/20);
figure
imagesc(F1>thr);part_on
title('ROI');

figure
imagesc(reshape(abs(kspace(40,:,:)),[Nx Ny]));
title('Slice 40 - K-space as directly acquired');
%kdat=kraw(1:2:end-1,:)+1i*kraw(2:2:end,:); idat=fftshift(ifft2(fftshift(kdat)));

figure
X = radon(abs(F1),0:180);
subplot(2,1,1);imagesc(iradon(X,0:180));title('Take Radon then invert');
subplot(2,1,2);imagesc(X);title('Sinogram');
%imagesc(X)

