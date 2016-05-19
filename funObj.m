function [loss, df ] = funObj( phase, intensity, z, nfocus, thresholdh, thresholdl, lambda, ps, Nx, Ny, maskFun, useGPU)
%FUNOBJ Summary of this function goes here
%   Detailed explanation goes here

phase = reshape(phase, [Nx, Ny]);
objectField = intensity * exp(1i * phase);
%objectField = phase;
loss = 0;

if useGPU
df = zeros(Nx, Ny, 'gpuArray');
objectField = gpuArray(objectField);
else
df = zeros(Nx, Ny);
end


for i = 1 : numel(z)
    
    HStack = GenerateFresnelPropagationStack(Nx, Ny, z(i)-z(nfocus), lambda, ps);
    mask = maskFun(z(i));    
    imagez = fftshift(fft2(objectField .* HStack));
    imageInten = abs(imagez.^2);
    maskh = mask .* (imageInten < thresholdh);
    maskl = (1-mask) .* (imageInten > thresholdl);
    
    diffh = maskh .* (imageInten - thresholdh);
    diffl = maskl .* (imageInten - thresholdl);
    
    temph = imagez.*diffh;
    temph = conj(HStack).*(Nx*Ny*ifft2(ifftshift(temph)));
    templ = imagez.*diffl;
    templ = conj(HStack).*(Nx*Ny*ifft2(ifftshift(templ)));
%     templ = Nx*Ny *abs(HStack.^2).* (objectField.*diffl);
%     temph = Nx*Ny *abs(HStack.^2).* (objectField.*diffh);

    loss = loss + sum(sum(diffh.^2 + diffl.^2)); 
    
    
    df = df +  temph + templ;
    %clear HStack mask imagez imageInten maskh maskl diffh diffl temph templ
end
%df = df .* (1i * intensity * exp(1i*phase));
df = - real(df).*sin(phase) + imag(df) .* cos(phase);
df = df(:);

loss = gather(loss);
df = gather(df);
end

