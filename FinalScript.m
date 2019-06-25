clear

SNR = @(image,truth) 10^(snr(abs(image),abs(image)-abs(truth))/20);
res = 256;

% Take fft(phantoms)
X1 = imgaussfilt(phantom); K1 = fftshift(fft2(phantom));
A1 = ifft2(K1);
figure;subplot(2,2,1);imagesc(abs(X1));title('truth');
subplot(2,2,2);imagesc(abs(A1));title('reconstructed');
subplot(2,2,3);imagesc(abs(X1)-abs(A1));title('noise');

SNR(A1,X1) %10 empirically

sigma = 1; %Experimentally, this gives a SNR of ~60
noise = normrnd(0,sigma,size(K1)) + 1i * normrnd(0,sigma,size(K1));

K2 = K1+noise;
A2 = ifft2(K2);
figure;subplot(2,2,1);imagesc(abs(X1));title('truth');
subplot(2,2,2);imagesc(abs(A1));title('reconstructed');
subplot(2,2,3);imagesc(abs(X1)-abs(A1));title('noise');

SNR(A2,A1) %10 empirically

%Create correlated noise
[xsp,ysp] = meshgrid(-128:127,-128:127);
quadrance = xsp.^2 + ysp.^2;
index = union(find(quadrance < 15^2),find(quadrance > 64^2));

N1 = unifrnd(-1,1,256,256);
N1(index) = 0;

N2real = imgaussfilt(real(ifft2(N1)),5);
N2imag = imgaussfilt(imag(ifft2(N1)),5);
N2 = N2real + 1i* N2imag;
%figure;imagesc(abs(IMAGE));


K3 = K1+1000000*N2;
A3 = ifft2(K3);
figure;subplot(2,2,1);imagesc(abs(X1));title('truth');colorbar;
subplot(2,2,2);imagesc(abs(A3));title('reconstructed');colorbar;
subplot(2,2,3);imagesc(abs(X1)-abs(A3));title('noise');

SNR(A3,A1) %10 empirically


% Apply Hermitian symmetry
K4 = zeros(size(K3));
K4(:,129:256) = K3(:,129:256);
K4(:,128:-1:1) = conj(K3(256:-1:1,129:256));
A4 = ifft2(K4);
figure;subplot(2,2,1);imagesc(abs(X1));title('truth');colorbar;
subplot(2,2,2);imagesc(abs(A4));title('reconstructed');colorbar;
subplot(2,2,3);imagesc(abs(X1)-abs(A4));title('noise');

SNR(A4,A1) %10 empirically


% Apply Hermitian symmetry
K5 = zeros(size(K2));
K5(:,129:256) = K2(:,129:256);
K5(:,128:-1:1) = conj(K2(256:-1:1,129:256));
A5 = ifft2(K5);
figure;subplot(2,2,1);imagesc(abs(X1));title('truth');colorbar;
subplot(2,2,2);imagesc(abs(A5));title('reconstructed');colorbar;
subplot(2,2,3);imagesc(abs(X1)-abs(A5));title('noise');

SNR(A5,A1) %10 empirically


kspace = fftshift(fft2(imgaussfilt(phantom)));

A6 = homodyneX(kspace,20,'step');
figure;subplot(2,2,1);imagesc(abs(X1));title('truth');colorbar;
subplot(2,2,2);imagesc(abs(A6));title('reconstructed - step');colorbar;
subplot(2,2,3);imagesc(abs(X1)-abs(A6));title('noise');


A7 = homodyneX(kspace,20,'ramp');
figure;subplot(2,2,1);imagesc(abs(X1));title('truth');colorbar;
subplot(2,2,2);imagesc(abs(A7));title('reconstructed - ramp');colorbar;
subplot(2,2,3);imagesc(abs(X1)-abs(A7));title('noise');


%Not undersampled
A8 = filteredBackProjection(X1,kspace,1.0,402);
A9 = filteredBackProjection(X1,kspace,0.5,402);

% -------- GRIDDING ---------
FOV = 1; %.28 in the code, but let's normalize.
kReach = res/2;

% Some preprocessing
bottom = -res/2; top = (res/2)-1;
[kXgrid,kYgrid] = meshgrid((bottom:top)/FOV,(bottom:top)/FOV);
tau = .15;
nearest = 8;


fsamp = .25; R = .25;
[T1, kres] = cartesianTraj(kReach,fsamp,R,FOV);
A10 = griddingTriangular(K1,T1,res,kres,nearest,FOV);
A11 = griddingGaussian(K1,T1,res,kres,kXgrid,kYgrid,tau,nearest,FOV);

disp('completed cartesian\n');

rSteps = 1.3*res;thSteps = 400;
[T2, kres] = polarTraj(kReach,rSteps,thSteps,FOV);
A12 = griddingTriangular(K1,T2,res,kres,nearest,FOV);
A13 = griddingGaussian(K1,T2,res,kres,kXgrid,kYgrid,tau,nearest,FOV);

disp('completed polar\n');

numSteps = 800;numLoops = 100;numLeaves = 100;
[T3, kres] = spiralTraj(kReach,numSteps,numLoops,numLeaves,FOV);
A14 = griddingTriangular(K1,T3,res,kres,nearest,FOV);
A15 = griddingGaussian(K1,T3,res,kres,kXgrid,kYgrid,tau,nearest,FOV);

figure;
subplot(3,3,1);imagesc(abs(A8));
subplot(3,3,2);imagesc(abs(A9));
subplot(3,3,3);imagesc(abs(A10));
subplot(3,3,4);imagesc(abs(A11));
subplot(3,3,5);imagesc(abs(A12));
subplot(3,3,6);imagesc(abs(A13));
subplot(3,3,7);imagesc(abs(A14));
subplot(3,3,8);imagesc(abs(A15));

% Apply noise to k-space.

% ifft to a baseline SNR of 60.

% exploit conjugate symmetry (??).
% Undersample (in one or both directions)
% Homodyne (take various strips in the middle, with and without noise)
% Filtered backprojection
% Gridding (polar/spiral)


% Do a wrapper for normal, homodyne, FBP, and gridding
% Do other phantoms
% Add noise in k-space. Calibrate it by putting an SNR of 50 on the
% "normal" one



% Look at scanned data
% IFFT(cartesian)
% HOMODYNE
% FBP