function [loss, df ] = SourceFunObj( phase_source, z, nfocus, lambda, ps, Nx, Ny, Ividmeas, useGPU, ratio)
%FUNOBJ Summary of this function goes here
%   Detailed explanation goes here

phase = reshape(phase_source(1:Nx*Ny), [Nx, Ny]);
source = reshape(phase_source(Nx*Ny+1:end), [Nx, Ny]);
objectField = exp(1i * phase);

%objectField = phase;
loss = 0;


if useGPU
dfsource = zeros(Nx, Ny, 'gpuArray');
dfphase = zeros(Nx, Ny, 'gpuArray');
objectField = gpuArray(objectField);
else
dfsource = zeros(Nx, Ny);
dfphase = zeros(Nx, Ny);
end


for i = 1 : numel(z)
    
    HStack = GenerateFresnelPropagationStack(Nx, Ny, z(i)-z(nfocus), lambda, ps);
    fieldz = fftshift(fft2(objectField .* HStack));
    coherent_spectral = fft2(ifftshift(abs(fieldz.^2)));
    
    source_spectral = fft2(source);
    imagez = ifft2(source_spectral .* coherent_spectral);
    
    diff = imagez - Ividmeas(:,:,i);
    temp_source = 2 * ifft2(conj(coherent_spectral) .* fft2(diff));
    temp_phase = fftshift(ifft2( conj(source_spectral).* fft2(diff)));
    temp_phase = fieldz .* temp_phase;
    temp_phase = conj(HStack).*(Nx*Ny*ifft2(ifftshift(temp_phase)));
%     templ = Nx*Ny *abs(HStack.^2).* (objectField.*diffl);
%     temph = Nx*Ny *abs(HStack.^2).* (objectField.*diffh);

    loss = loss + sum(sum(diff.^2)); 
    
    dfphase = dfphase + temp_phase;
    dfsource = dfsource + temp_source;
    %clear HStack mask imagez imageInten maskh maskl diffh diffl temph templ
end
%df = df .* (1i * intensity * exp(1i*phase));

df = zeros(2*Nx*Ny, 1);
loss = gather(loss);
dfphase = -real(dfphase) .* sin(phase) + imag(dfphase).*cos(phase);
df(1:Nx*Ny) = gather(dfphase(:)) * ratio;
%df(Nx*Ny+1:end) = gather(dfsource(:)) / ratio;

end

