function LV_competition_function_jitter_grid(L,N,D,P,pillarq,R,dx,rep,jitt_R,jitt_dx,jitt_num,seed)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nick Lowery and Tristan Ursell
% 2018
%
% Structured environments fundamentally alter dynamics and stability of ecological communities
% https://www.biorxiv.org/content/early/2018/07/10/366559
% 
% This function runs a non-dimensionalized modified 2D Lotka-Volterra ecological
% simulation of three species competing in an intransitive loop, and outputs 
% images of the simulation and a file containing the associated metadata. It 
% requires the image_save.m script to be in the active directory. For these 
% simulations, the pillar lattice is provided pre-made using the pillar_jitter.m
% function, which adds random noise to either the pillar sizes and/or spacing.
% These grid (.mat) files should reside in the current directory.
%
% Input parameters:
% L = size of the simulation box (pixels)
% N = number of simulation iterations (dt units, specified in the script)
% D = diffusion coefficient (dimensionless)
% P = average interaction strength (dimensionless)
% pillarq = logical (1,0); should pillars be included?
% R = pillar radius (pixels)
% dx_mult = center-center pillar spacing (units of R)
% rep = unique replicate identifier
% jitt_R, jitt_dx, jitt_num = identifiers for specific pillar lattice file
% seed = random seed for density matrix initialization
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%*******************************************
%********* INITIALIZATION ******************
%*******************************************

rng(seed)

%output images?
imageq = 1;

%time step in ND units
dt = 0.01;

%local concentration discrete cutoff (will be set to 0 below this value)
dsc = 0.01;

%frequency to show results
showt = 100;

%intialize the density matrices
A = 1/10 * rand(L,L);
B = 1/10 * rand(L,L);
C = 1/10 * rand(L,L);

%diffusion convolution filter
sigma = sqrt(4*D*dt);
if mod(ceil(3*sigma),2) == 0
    win_sz = ceil(3*sigma) + 1;
else
    win_sz = ceil(3*sigma);
end
Gauss = fspecial('gaussian', [win_sz,win_sz], sigma);
rescale_filt = 1 ./ conv2(ones(size(A)), Gauss, 'same');

if L < win_sz
    warning('Diffusion window size is larger than simulation box.')
end

%file IO
if imageq == 1
    path0 = cd;
    
    file1 = ['3sp_LV_comp_jitter_grid_L-' num2str(L) ...
	'_N-' num2str(N) ...
        '_D-' num2str(D) ...
        '_P-' num2str(P) ...
	'_pillars-' num2str(pillarq) ...
        '_R-' num2str(R) ...
        '_dx-' num2str(dx) ...
        '_dsc-' num2str(dsc) ...
        '_rep-' num2str(rep) ...
	'_jittR-' num2str(jitt_R) ...
	'_jittdx-' num2str(jitt_dx) ...
	'_jittnum-' num2str(jitt_num)];
    
    path1 = fullfile(path0,file1);
    mkdir(path1)
    
    q = 0;
end

%tracking vectors (global mean, scale factor for each image)
meanA = zeros(N,1);
meanB = zeros(N,1);
meanC = zeros(N,1);
scale_fac = zeros(N,3);

%*******************************************
%********* PILLARS *************************
%*******************************************

% read in pre-defined jittered grid
    % pre-screened for sufficient inter-pillar distance for diffusion to work normally
grids = dir('jitter*.mat');
grid_names = {grids(:).name};
str = ['R-' num2str(R) '_dx-' num2str(dx/R) '_vR-' num2str(jitt_R) '_vdx-' num2str(jitt_dx) '_num-' num2str(jitt_num)];
filt_all_temp = load(string(grid_names(find(contains(grid_names, str)))));
filt_all = double(filt_all_temp.filt_all);
   
%sum weights of active areas for calcualting mean abundances
filt_all_weight = sum(filt_all(:));
    
%prepare reflecting boundary conditions in grid
conv_filt = conv2(filt_all,Gauss,'same');
rescale_filt = filt_all ./ conv_filt;
rescale_filt(isnan(rescale_filt)) = 0;

%*******************************************
%********** SIMULATION *********************
%*******************************************

for i = 1:N
    %diffusion convolutions
    A = conv2(A, Gauss, 'same') .* rescale_filt;
    B = conv2(B, Gauss, 'same') .* rescale_filt;
    C = conv2(C, Gauss, 'same') .* rescale_filt;
        
    %interspecies interactions    
    dA = dt * A .* (1 - (A+B+C)) .* (1 - C/P);
    dB = dt * B .* (1 - (A+B+C)) .* (1 - A/P);
    dC = dt * C .* (1 - (A+B+C)) .* (1 - B/P);
    
    %deterministic update populations
    A = A + dA;
    B = B + dB;
    C = C + dC;
    
    %hard upper bound at carrying ccapacity
    A(A > 1) = 1;
    B(B > 1) = 1;
    C(C > 1) = 1;
    
    %discrete cutoff (hard lower bound)
    A(A < dsc) = 0;
    B(B < dsc) = 0;
    C(C < dsc) = 0;
    
    %calculate tracking vectors
    if pillarq == 1
        meanA(i) = sum(A(filt_all == 1)) / filt_all_weight;
        meanB(i) = sum(B(filt_all == 1)) / filt_all_weight;
        meanC(i) = sum(C(filt_all == 1)) / filt_all_weight;
    else
        meanA(i) = mean(A(:));
        meanB(i) = mean(B(:));
        meanC(i) = mean(C(:));
    end
    
    %scale factor - maximize contrast in output images
    	%(use to revert to true value for later quantitiative analysis)
    scale_fac(i,1) = max(A(:));
    scale_fac(i,2) = max(B(:));
    scale_fac(i,3) = max(C(:));
    
    %output images
    if mod(i,showt) == 1
        disp(['step ' num2str(i) ' of ' num2str(N) ' for ' file1])
	drawnow('update')        
 
        %grey filter
        grey_filt = (A == 0) .* (B == 0) .* (C == 0);
        
        Atemp = A; 
        Btemp = B;
        Ctemp = C;
           
        if pillarq == 1
            Atemp(grey_filt == 1) = scale_fac(i,1) / 2;
            Btemp(grey_filt == 1) = scale_fac(i,2) / 2;
            Ctemp(grey_filt == 1) = scale_fac(i,3) / 2;
        end
        
        im_out(:,:,1) = uint8(Atemp / scale_fac(i,1) * 255);
        im_out(:,:,2) = uint8(Btemp / scale_fac(i,2) * 255);
        im_out(:,:,3) = uint8(Ctemp / scale_fac(i,3) * 255);
        
        if imageq == 1
            q = q + 1;
            image_save(im_out, fullfile(path1, [file1 '_' sprintf('%04d',q) '.png']))
        end
    end
end

%*******************************************
%********** METADATA ***********************
%*******************************************

meanA_out = meanA(1:showt:N);
meanB_out = meanB(1:showt:N);
meanC_out = meanC(1:showt:N);

scale_fac_out = scale_fac(1:showt:N, :)

save([file1 '_metadata.mat'], 'L', 'N', 'D', 'P', 'pillarq', 'R', 'dx', 'dt', 'dsc', ...
	'meanA_out', 'meanB_out', 'meanC_out', 'scale_fac_out', 'rep', 'seed', ...
	'jitt_R', 'jitt_dx', 'jitt_num', 'filt_all_weight')

disp(['Done with ' file1 '.'])

%exit

end
