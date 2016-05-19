%% Setup params
% All length value has unit meter in this file.
% The 3d region is behind lens after SLM. 
% Since there is tube lens and objective lens after the the region. We make
% the resolution 20 times that of the region after objective lens to
% account for the demagnification.
addpath(genpath('minFunc_2012'))

resolutionScale = 20; % The demagnification scale of tubelens and objective.
lambda = 1e-6;  % Wavelength
psHolograph = 2e-6 * resolutionScale;      % Pixel Size (resolution) at the scattered 3D region
psSLM = 40e-6;      % Pixel Size (resolution) at the scattered 3D region
Nx = 300;       % Number of pixels in X direction
Ny = 400;       % Number of pixels in Y direction
useGPU = 0;     % Use GPU to accelerate computation. Effective when Nx, Ny is large (e.g. 600*800).
focal_SLM = 0.2; % focal length of the lens after slm.
% Please note that this program does not take into account of subsampling
% issues. Instead, it assumes that the SLM has the same resolution and 
% size as the reconstructed images in 3D region. i.e. the SLM is right
% before the objective lens. I think this assumption can be satisfied by
% magnify/minify SLM using an optical system. We can discuss if this needs
% to be changed. (Jingzhao)

z = [400 :4: 600] * 1e-6 * resolutionScale;   % Depth level requested in 3D region.
nfocus = 20;                % z(nfocus) denotes the depth of the focal plane.
thresholdh = 20000000;          % Intensity required to activate neuron.
thresholdl = 0;             % Intensity required to not to activate neuron.

%% Point Targets
radius = 10 * 1e-6 * resolutionScale; % Radius around the point.
targets = [0,0,450; 100, 100, 500; -100,-200,550;] * 1e-6 * resolutionScale; % Points where we want the intensity to be high.
maskfun = @(zi)  generatePointMask( targets, radius, zi, Nx, Ny, psHolograph, useGPU);


%% Complex Target
% load('complexAB');
% zrange1 = [450,480] * 1e-6;
% zrange2 = [550,580] * 1e-6;
% 
% maskfun = @(zi) generateComplexMask( zi, Nx, Ny, maskA, zrange1, maskB, zrange2);



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
tic;


% Scale the gradient of phase by 1 and the gradient of source by 0.
% This makes sure that only the phase is updated in each iteration.
ratio_phase = 1;
ratio_source = 0;


f = @(x)SourceFunObj(x, z, nfocus, lambda, psSLM, focal_SLM, Nx, Ny, thresholdh, thresholdl, maskfun,  useGPU, ratio_phase, ratio_source);
%phase_source = minFunc(f, x0, options);


matlab_options = optimoptions('fmincon','GradObj','on', 'display', 'iter', ...
    'algorithm','interior-point','Hessian','lbfgs', 'MaxIter', 30);
lb = -inf(2*Nx*Ny, 1);
lb(end/2+1:end) = 0;
ub = inf(2*Nx*Ny, 1);
nonlcon = [];
phase_source = fmincon(f,x0,[],[],[],[],lb,ub,nonlcon,matlab_options);


% The following part optimizes phase and source at the same time.
% ratio_phase = 1;
% ratio_source = 1; 
% f = @(x)SourceFunObj(x, z, nfocus, lambda, psSLM, focal_SLM, Nx, Ny, thresholdh, thresholdl, maskfun,  useGPU, ratio_phase, ratio_source);
% %phase_source = minFunc(f, phase_source, options);
% phase_source = fmincon(f,phase_source,[],[],[],[],lb,ub,nonlcon,matlab_options);
toc;

phase = reshape(phase_source(1:Nx*Ny), [Nx, Ny]);
source = reshape(phase_source(Nx*Ny+1:end), [Nx, Ny]);

%% plot
 
Ividmeas = zeros(Nx, Ny, numel(z));
figure();
for i = 1:numel(z)
    HStack = GenerateFresnelPropagationStack(Nx,Ny,z(i) - z(nfocus), lambda, psSLM, focal_SLM);
    imagez = fresnelProp(phase, source, HStack);
    Ividmeas(:,:,i) = imagez;
    imagesc(imagez);colormap gray;colorbar;title(sprintf('Distance z %d', z(i)));
    caxis([0, 5e6]);
    pause(0.1);
end
