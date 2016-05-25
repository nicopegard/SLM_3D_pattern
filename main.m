clear all;
tag = 'slmFocal_20_600by800_AB';


%% Setup params
% All length value has unit meter in this file.
% The 3d region is behind lens after SLM. 

addpath(genpath('minFunc_2012'))

resolutionScale = 20; % The demagnification scale of tubelens and objective. f_tube/f_objective
lambda = 1e-6;  % Wavelength
focal_SLM = 0.2; % focal length of the lens after slm.
psSLM = 20e-6;      % Pixel Size (resolution) at the scattered 3D region
Nx = 600;       % Number of pixels in X direction
Ny = 800;       % Number of pixels in Y direction


psXHolograph = lambda * focal_SLM/ psSLM / resolutionScale / Nx;      % Pixel Size (resolution) at the scattered 3D region
psYHolograph = lambda * focal_SLM/ psSLM / resolutionScale / Ny;      % Pixel Size (resolution) at the scattered 3D region

useGPU = 1;     % Use GPU to accelerate computation. Effective when Nx, Ny is large (e.g. 600*800).


z = [400 :4: 600] * 1e-6 ;   % Depth level requested in 3D region.
nfocus = 20;                % z(nfocus) denotes the depth of the focal plane.
thresholdh = 20000000;          % Intensity required to activate neuron.
thresholdl = 0;             % Intensity required to not to activate neuron.

%% Point Targets
radius = 10 * 1e-6 ; % Radius around the point.
targets = [150,150,450; 0, 0, 500; -150,-150,550;] * 1e-6 ; % Points where we want the intensity to be high.
maskfun = @(zi)  generatePointMask( targets, radius, zi, Nx, Ny, psXHolograph,psYHolograph, useGPU);


%% Complex Target
% load('largeAB');
% zrange1 = [450,480] * 1e-6;
% zrange2 = [550,580] * 1e-6;
% 
% maskfun = @(zi) generateComplexMask( zi, Nx, Ny, maskA, zrange1, maskB, zrange2);
% 


%% Optimization

% options.Method = 'lbfgs';
% options.optTol = 10^(-6);
% options.progTol = 10^(-6);
% options.MaxIter = 30;
% options.MaxFunEvals = 350;


% The starting point. reshape(x0(1:Nx*Ny), [Nx, Ny]) encodes the phase on
% SLM in rads. Normally initialized to zeros. reshape(x0(1+Nx*Ny:end), [Nx, Ny])
% encodes the source intensity. Need to be nonnegative.
x0 = ones(2*Nx*Ny, 1) * 1e-10;
% This sets a coherent light source.
x0(end/2 + Nx*(Ny/2 +0.5)) = 1;
% x0(Nx*Ny+1:end) = randn([Nx*Ny, 1])/Nx + 1/Nx/Ny;
% x0 = x0 .* (x0>0) + 1/Nx/Ny *(x0<0);

tic;


% Scale the gradient of phase by 1 and the gradient of source by 0.
% This makes sure that only the phase is updated in each iteration.
ratio_phase = 1;
ratio_source = 0;


HStacks = zeros(Nx, Ny, numel(z));
if useGPU
    HStacks = gpuArray(HStacks);
end
for i = 1 : numel(z)
    HStacks(:,:,i) = GenerateFresnelPropagationStack(Nx, Ny, z(i)-z(nfocus), lambda, psXHolograph,psYHolograph, useGPU);
end
kernelfun = @(x) HStacks(:,:,x);
% kernelfun = @(i) GenerateFresnelPropagationStack(Nx, Ny, z(i)-z(nfocus), lambda,psXHolograph,psYHolograph, focal_SLM, useGPU);


f = @(x)SourceFunObj(x, z, Nx, Ny, thresholdh, thresholdl, maskfun, kernelfun, useGPU, ratio_phase, ratio_source);
%phase_source = minFunc(f, x0, options);


matlab_options = optimoptions('fmincon','GradObj','on', 'display', 'iter', ...
    'algorithm','interior-point','Hessian','lbfgs', 'MaxFunEvals', 50, 'MaxIter', 30);
lb = -inf(2*Nx*Ny, 1);
lb(end/2+1:end) = 0;
ub = inf(2*Nx*Ny, 1);
nonlcon = [];
phase_source1 = fmincon(f,x0,[],[],[],[],lb,ub,nonlcon,matlab_options);


phase1 = reshape(phase_source1(1:Nx*Ny), [Nx, Ny]);
source1 = reshape(phase_source1(Nx*Ny+1:end), [Nx, Ny]);
hologram = floor(mod(phase1, 2*pi)/2/pi * 255);

f = figure(1)
imagesc(hologram);
title('SLM hologram');

%% The following part optimizes phase and source at the same time.
% ratio_phase = 1;
% ratio_source = 1; 
% f = @(x)SourceFunObj(x, z, Nx, Ny, thresholdh, thresholdl, maskfun, kernelfun, useGPU, ratio_phase, ratio_source);
% %phase_source = minFunc(f, phase_source, options);
% phase_source2 = fmincon(f,phase_source1,[],[],[],[],lb,ub,nonlcon,matlab_options);
% 
% 
% phase2 = reshape(phase_source2(1:Nx*Ny), [Nx, Ny]);
% source2 = reshape(phase_source2(Nx*Ny+1:end), [Nx, Ny]);

toc;

%% plot
Ividmeas = zeros(Nx, Ny, numel(z));
usenoGPU = 0;
figure();
for i = 1:numel(z)
    HStack = GenerateFresnelPropagationStack(Nx,Ny,z(i) - z(nfocus), lambda, psXHolograph,psYHolograph, usenoGPU);
    imagez = fresnelProp(phase1, source1, HStack);
    Ividmeas(:,:,i) = imagez;
    imagesc(imagez);colormap gray;title(sprintf('Distance z %d', z(i)));
    colorbar;
    caxis([0, 5e6]);
%     filename = sprintf('pointTarget%d.png', i);
%     print(['data/' filename], '-dpng')
    pause(0.1);
end
%save(['source_phase_result_' tag '.mat'], 'source1', 'phase1', 'source2', 'phase2', 'hologram');





