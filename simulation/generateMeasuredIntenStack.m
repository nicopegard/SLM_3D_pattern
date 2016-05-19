%% save
z = [200 :16: 1000] * 1e-6;
nfocus = 20;

%phase = 2*randi(300,300) - 1;

measuredInten = zeros(Nx, Ny, numel(z));
figure();
%objectField = exp(1i*2*phase)/1000;
TrueImg = TrueImg(1:300, 1:300);
objectField = TrueImg(1:300, 1:300);
for i = 1:numel(z)
    HStack = GenerateFresnelPropagationStack(Nx, Ny, z(i)-z(nfocus), lambda, ps, useGPU);
    imagez = fftshift(fft2(objectField .* HStack));
    imageInten = abs(imagez.^2);
    measuredInten(:,:,i) = imageInten;
    imagesc(imageInten);colormap gray;colorbar;title(sprintf('Distance z %d', z(i)));
    %caxis([0, 1e2]);
    pause(0.1)


end
truephase = phase;

save('StarMeasuredintenStack.mat','z', 'nfocus','lambda', 'ps', 'Nx', 'Ny','measuredInten', 'TrueImg');