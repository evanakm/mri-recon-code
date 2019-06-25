%clear
%filename = strcat('spgr3d_256x_256y_160z_1e_12c-slice',num2str(slice),'.fft');

%fileID = fopen(filename);
%XX = fread(fileID,'uint16');

%mat = fftshift(reshape(X,[256 256]));
%mat = reshape(XX,[256 256]);
mat = fftshift(phantom);

%fclose(fileID);

figure;
truth=abs(fftshift(mat));
imagesc(truth);
title('Reconstructed image');

noise = normrnd(0,.50,size(mat)); %Experimentally, this gives a SNR of ~45
SNR = 10^(snr(mat,noise)/20);

thr = 25;
roimat = mat(mat>thr); roinoise = noise(mat>thr);

ROIsnr = 10^(snr(roimat,roinoise)/20)

%figure;
%imagesc(mat>thr)

%add noise to truth.
mat = mat+noise;

%kspace = fft2(mat)/(size(mat,1)*size(mat,2));
kspace = fft2(mat);
%kspace = fft2(truth)/(size(truth,1)*size(truth,2));
%kspace = fft2(fftshift(truth))/(size(truth,1)*size(truth,2));
KSP = fftshift(kspace);
%KSP = kspace;
mag = abs(KSP);
%figure;imagesc(mag);colorbar;title('k-space(cartesian)')

angles = 402;
%polarK = zeros(256,402);
thDEG = (0:(angles-1))*180/angles;
thRAD = (0:(angles-1))*pi/angles;

SIZE = ceil(256*sqrt(2));
if(mod(SIZE,2)==0)
    SIZE = SIZE+1; %Make sure SIZE is odd
end
BASE = (SIZE-1)/2;

%SIZE = 512;BASE = 256;

polarK = zeros(SIZE,angles);

for j = 1:angles
    for i = 1:SIZE
%       r = (-BASE-1+i)/sqrt(2);
       r = (-BASE-1+i);
       th = thRAD(j);

       x = 129+r*cos(th);
       y = 129+r*sin(th);
%       x = 257+r*cos(th);
%       y = 257+r*sin(th);
       
       if(x>=1 && x<257 && y>=1 && y<257)
           polarK(i,j) = squareinterp(x,y,KSP);       
       else
           polarK(i,j) = 0;
       end
       
    end
end


%for j = 1:angles
%    for i = 1:256
%       r = (-BASE-1+i)/sqrt(2);
%       r = (-BASE-1+i);
%       th = thRAD(j);
%       x = 129+r*cos(th);
%       y = 129+r*sin(th);
%                   
%       polarK(i,j) = squareinterp(x,y,KSP);
%       
%    end
%end

%figure;imagesc(abs(polarK));colorbar;title('k-space(polar)')


%idat=fftshift(ifft(fftshift(polarK,1),[],1),1);
%idat=fftshift(ifft(polarK,[],1),1);
%idat=ifft(polarK,[],1);
%idat=fftshift(ifft(polarK,[],1),2);%best one
%idat=fftshift(ifft(polarK,[],1));%out of phase but gets details
%idat=fftshift(ifft(polarK,[],1));

idat=fftshift(ifft(polarK,[],1),1);
%idat=ifft(polarK,[],1);
idat=[idat(2:end,:); idat(1,:)];
%m=iradon(abs(idat),thDEG,'linear','Ram-Lak',1,512);
m=iradon(abs(idat),thDEG,'linear','ram-lak',1.0,floor(256*sqrt(2)));
%m=iradon(abs(idat),thDEG);

SZ = size(m,2);
finalIM = abs(m');
finalIM = finalIM(:,SZ:-1:1);
finalIM = finalIM/max(finalIM(:)); %Normalize
finalIM = imresize(finalIM,[256 256]);

%finalIM = finalIM .*(truth>0); %Remove outside ROI

sgrm = radon(finalIM,thDEG);
%middle = ceil(size(sgrm,1)/2);
%range = (middle-183):(middle+183);
%sgrm = sgrm(range,:);

figure
subplot(2,1,1)
imagesc(finalIM);title('From iradon');
subplot(2,1,2)
%imagesc(abs(idat));title('Sinogram')
imagesc(sgrm);title('Sinogram')


Z1 = sgrm/max(sgrm(:));
%int = X(3:365,:);
Z2 = abs(X)/max(X(:));

RADIAL = size(Z1,1);
Z2 = imresize(Z2,[RADIAL angles]);

figure
subplot(2,2,1)
imagesc(Z2);colorbar;title('Direct Radon');
subplot(2,2,2)
imagesc(Z1);colorbar;title('Reconstructed Radon');
subplot(2,2,3)
imagesc(Z1-Z2);colorbar;title('Difference');
subplot(2,2,4);
imagesc(iradon(Z1-Z2,thDEG));title('Image Difference');




%finalSNR = 10^(snr(truth,finalIM-truth)/20)

%figure
%imagesc(mat');

