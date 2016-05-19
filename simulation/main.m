%% Setup params
% All length value has unit meter in this file.
addpath(genpath('minFunc_2012'))

lambda = 1e-6;  % Wavelength
ps = 2e-6;      % Pixel Size (resolution) at the scattered 3D region
Nx = 300;       % Number of pixels in X direction
Ny = 300;       % Number of pixels in Y direction
useGPU = 0;     % Use GPU to accelerate computation. Effective when Nx, Ny is large (e.g. 600*800).

% Please note that this program does not take into account of subsampling
% issues. Instead, it assumes that the SLM has the same resolution and 
% size as the reconstructed images in 3D region. i.e. the SLM is right
% before the objective lens. I think this assumption can be satisfied by
% magnify/minify SLM using an optical system. We can discuss if this needs
% to be changed. (Jingzhao)


%% Point Targets
z = [400 :4: 600] * 1e-6;   % Depth level requested in 3D region.
nfocus = 20;                % z(nfocus) denotes the depth of the focal plane.

% simulation init
load('starTarget.mat');
load('complextarget.mat');
phase = mask;
source = star(1:Nx, 1:Ny);

%% Optimization

%load('SourceRecoverySimulationData.mat');
load('SourceRecoveryPoints');

% options = struct('GradObj','on','Display','iter','HessUpdate','bfgs','InitialHessType','identity','GoalsExactAchieve',0,...
%     'TolFun', 1e-12, 'TolX', 1e-12);
% [x2,fval2,exitflag,output,grad] = fminlbfgs(f,phase,options);
options.Method = 'lbfgs';
options.optTol = 10^(-6);
options.progTol = 10^(-6);
options.MaxIter = 100;
options.MaxFunEvals = 350;
x0 = zeros(2*Nx*Ny, 1);
x0(end * 3/4) = 1;
%x0(1:Nx*Ny) = phase(:);
%x0(1+Nx*Ny:end) = source(:);
tic;
f = @(x)SourceFunObj( x, z, nfocus, lambda, ps, Nx, Ny, Ividmeas, useGPU, 1);
phase_source = minFunc(f, x0, options);
% f = @(x)SourceFunObj( x, z, nfocus, lambda, ps, Nx, Ny, Ividmeas, useGPU, 1);
% phase_source = minFunc(f, phase_source, options);

toc;
phase = reshape(phase_source(1:Nx*Ny), [Nx, Ny]);
source = reshape(phase_source(Nx*Ny+1:end), [Nx, Ny]);







%% plot

% simulation init
load('pointHologram.mat');
load('complextarget.mat');
phase = phase;
source = zeros(Nx, Ny);
source(151,151) = 1;

Ividmeas = zeros(Nx, Ny, numel(z));
figure();
for i = 1:numel(z)
    HStack = GenerateFresnelPropagationStack(Nx,Ny,z(i) - z(nfocus), lambda, ps);
    imagez = fresnelProp(phase, source, HStack);
    Ividmeas(:,:,i) = imagez;
    imagesc(imagez);colormap gray;colorbar;title(sprintf('Distance z %d', z(i)));
    %caxis([0, 4e2]);
    pause(0.1);

end
save('SourceRecoveryPoints.mat');