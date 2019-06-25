kspace = fftshift(fft2(imgaussfilt(phantom)));

subplot(2,2,1);imagesc(homodyneX(kspace,5,'step'));
subplot(2,2,2);imagesc(homodyneX(kspace,5,'ramp'));
subplot(2,2,3);imagesc(homodyneX(kspace,5,'ramp')-homodyneX(kspace,5,'step'));





