function finalIM = filteredBackProjection(truth,kspace,options)

rStepSize = options.rStepSize;
angles = options.angles;
filter = options.filter;

[Nx,Ny] = size(kspace);

[X,Y] = meshgrid((-Nx/2):((Nx/2)-1),(-Nx/2):((Nx/2)-1));
Z = sqrt(X.^2+Y.^2);
Z = 0.5*max(Z(:)) - Z;
filt = max(Z,0);

%inds = 1:undersampleFactor:255;

%sliceUnder = zeros(floor(Nx/undersampleFactor,Ny/undersampleFactor);

%kdat1 = kdat + normrnd(0,1,[Nx,Ny])+1i*normrnd(0,1,[Nx,Ny]);
%kdat2 = kdat .* filt;
%kdat3 = kdat1 .* filt;

X = radon(abs(truth),0:180);
%figure;subplot(2,1,1);imagesc(iradon(X,0:180));title('Take Radon then invert');
%subplot(2,1,2);imagesc(X);title('Sinogram');

mat = fftshift(truth);

%kspaceRad = fft2(mat);KSP = fftshift(kspaceRad);
KSP = kspace;
thDEG = (0:(angles-1))*180/angles;
thRAD = (0:(angles-1))*pi/angles;

SIZE = ceil(Nx*sqrt(2));
if(mod(SIZE,2)==0)
    SIZE = SIZE+1; %Make sure SIZE is odd
end
BASE = (SIZE-1)/2;

R = (-BASE):rStepSize:BASE;
polarK = zeros(length(R),angles);

indices = zeros(length(R)*length(thRAD),2); %Debugging only
counter = 0;

for j = 1:angles
    for i = 1:length(R)

        r  = R(i);
        th = thRAD(j);

%        idxKX = trajectory(i,1)*FOV + Nx/2;
%        idxKY = trajectory(i,2)*FOV + Ny/2;

        idxKX = 1 + r*cos(th) + Nx/2;
        idxKY = 1 + r*sin(th) + Ny/2;
        
        counter = counter+1;
        indices(counter,:) = [idxKX, idxKY];
        


        fracX = idxKX - floor(idxKX);
        fracY = idxKY - floor(idxKY);

        
        if(idxKX<-0.5 || idxKX>=(Nx+0.5) || idxKY<-0.5 || idxKY>=(Ny+0.5))
            polarK(i,j) = 0;
        else
            
            % Convolving with a triangle
            if(fracX >= 0.5)
                lambdaX = 2 * fracX - 1; % reduces from 1 - (1-fracX)/.5, but this is not clear at first glance.
                xPT = ceil(idxKX);
            else
                lambdaX = 1 - 2 * fracX;
                xPT = floor(idxKX);
            end

            if(fracY >= 0.5)
                lambdaY = 2 * fracY - 1; % reduces from 1 - (1-fracX)/.5, but this is not clear at first glance.
                yPT = ceil(idxKY);
            else
                lambdaY = 1 - 2 * fracY;
                yPT = floor(idxKY);
            end

            if(xPT == 0)
                xPT = Nx;
            end

            if(yPT == 0)
                yPT = Ny;
            end

            polarK(i,j) = lambdaX*lambdaY*KSP(xPT,yPT);
            
        end
       
    end
end

idat=fftshift(ifft(polarK,[],1),1);
%idat=ifft(polarK,[],1);
idat=[idat(2:end,:); idat(1,:)];
%m=iradon(abs(idat),thDEG,'linear','Ram-Lak',1,512);
m=iradon(abs(idat),thDEG,'linear',filter,1.0,floor(Nx*sqrt(2)));
%m=iradon(abs(idat),thDEG);

SZ = size(m,2);
finalIM = fftshift(abs(m'));
finalIM = finalIM(:,SZ:-1:1);
finalIM = imresize(finalIM,[Nx Ny]);

% Originally returned normalized array, but some further analysis requires the
% unnormalized result.
normalizedIM = finalIM/max(finalIM(:)); %Normalize

sgrm = radon(normalizedIM,thDEG);

%figure
%subplot(2,1,1)
%imagesc(finalIM);title('From iradon');
%subplot(2,1,2)
%imagesc(abs(idat));title('Sinogram')
%imagesc(sgrm);title('Sinogram')


Z1 = sgrm/max(sgrm(:));
Z2 = abs(X)/max(X(:));

RADIAL = size(Z1,1);
Z2 = imresize(Z2,[RADIAL angles]);

if(options.comparePlots)
figure;title(strcat('FBP: rStepSize = ',num2str(rStepSize),', angles = ',num2str(angles)));
subplot(2,2,1);imagesc(Z2);colorbar;title('Direct Radon');
subplot(2,2,2);imagesc(Z1);colorbar;title('Reconstructed Radon');
subplot(2,2,3);imagesc(Z1-Z2);colorbar;title('Difference');
subplot(2,2,4);imagesc(iradon(Z1-Z2,thDEG));title('Image Difference');
end


%finalSNR = 10^(snr(truth,finalIM-truth)/20)

%figure
%imagesc(mat');


