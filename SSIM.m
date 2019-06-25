function [ whole,onlyROI ] = SSIM( image, reference, mask )
%SSIM Summary of this function goes here
%   Detailed explanation goes here

    whole = ssim(abs(image),abs(reference));
    
    image(mask) = 0;
    reference(mask) = 0;
    
    onlyROI = ssim(abs(image),abs(reference));
    


end

