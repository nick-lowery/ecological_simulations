function LV_competition_function(L,N,D,P,pillarq,R,dx_mult,rep,seed)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nick Lowery and Tristan Ursell
% 2018
%
% Structured environments fundamentally alter dynamics and stability of ecological communities
% https://www.biorxiv.org/content/early/2018/07/10/366559
% 
% This function runs a non-dimensionalized modified 2D Lotka-Volterra ecological
% simulation of a three species community competing in an intransitive loop, 
% and outputs images of the simulation and a file containing the associated
% metadata. It requires the image_save.m script to be in the active directory.
%
% Input parameters:
% L = size of the simulation box (pixels)
% N = number of simulation iterations (dt units, specified in the script)
% D = diffusion coefficient (dimensionless)
% P = interaction strength (dimensionless)
% pillarq = logical (1,0); should pillars be included?
% R = pillar radius (pixels)
% dx_mult = center-to-center pillar spacing (units of R)
% rep = UNIQUE replicate identifier
% seed = set seed for density matrix initialization
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
    
    file1 = ['3sp_LV_comp_L-' num2str(L) ...
	'_N-' num2str(N) ...
        '_D-' num2str(D) ...
        '_P-' num2str(P) ...
        '_pillars-' num2str(pillarq) ...
        '_R-' num2str(R) ...
        '_dx-' num2str(dx_mult*R) ...
        '_dsc-' num2str(dsc) ...
        '_rep-' num2str(rep)];
    
    path1 = fullfile(path0, file1);
    mkdir(path1)
    
    %image index tracking
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
if pillarq == 1
    %create spatial matrices / filters
    Xmat = ones(L,1)*(1:L);
    Ymat = Xmat';
    
    %create centers for obstruction pattern
    dx = dx_mult*R;
    
    %form hexagonal grid position arrays
    x_cent0 = 0:dx:(L+dx);
    y_cent0 = [0:(sqrt(3)/2*dx):(L+dx*sqrt(3)/2)] + R;
    
    x_cent = zeros(length(y_cent0), length(x_cent0));
    for i = 1:length(y_cent0)
        if mod(i,2) == 1
            x_cent(i,:) = x_cent0 + dx/4;
        else
            x_cent(i,:) = x_cent0 - dx/4;
        end
    end
    
    y_cent = zeros(size(x_cent));
    for j = 1:length(x_cent0)
        y_cent(:,j) = y_cent0';
    end
    
    %create pillars
    filt_all = ones(size(Xmat));
    for i = 1:length(x_cent0)
        for j = 1:length(y_cent0)
            filt_temp = ((Xmat - x_cent(j,i)) .^ 2 + (Ymat - y_cent(j,i)) .^ 2) >= R^2;
            filt_all = filt_all .* filt_temp;
        end
    end
    
    dbound = -L/2;
    filt_all_temp = zeros(size(filt_all));
    filt_all_temp(1:round(L/2) - dbound,:) = filt_all(1:round(L/2) - dbound,:);
    filt_all_temp(round(L/2) - dbound + 1:end,:) = 1;
    filt_all = filt_all_temp;
    
    %calculate weight of non-pillar area for calcualting mean abundances
    filt_all_weight = sum(filt_all(:));
    
    %prepare reflecting boundary conditions in grid
    conv_filt = conv2(filt_all, Gauss, 'same');
    rescale_filt = filt_all ./ conv_filt;
    rescale_filt(isnan(rescale_filt)) = 0;
else
    R = 0; dx = 0;
    filt_all = ones(L,L);
end

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

save([file1 '_metadata.mat'], 'L', 'N', 'D', 'P', 'pillarq', 'R', 'dx', 'dt', 'dsc', 'seed', ...
	'meanA_out', 'meanB_out', 'meanC_out', 'scale_fac', 'rep', 'filt_all_weight')

disp(['Done with ' file1 '.'])

%exit

end
