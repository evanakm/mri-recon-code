function snr = SNR(image,truth,ROI)

    image = abs(image);
    truth = abs(truth);
    noise = image - truth;
    
    [NX,NY] = size(image);
    
    truthInt = truth.^2;
    signal = max(truthInt(ROI));
    
    [indsX,indsY] = meshgrid(1:NX,1:NY);
    X = reshape(indsX,[numel(indsX),1]);
    Y = reshape(indsY,[numel(indsY),1]);
    
    nvec = reshape(noise,[numel(noise),1]);
    
    data = [ ones(size(X)) X Y ];
    [a,b,resid] = regress(nvec,data);
    
    noise = std(resid);
    
    snr = [];
    snr.snr = signal/noise;
    snr.signal = signal;
    snr.noise = noise;
    
    
    
end