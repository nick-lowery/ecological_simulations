function LV_competition_two_species_asym(L,N,D,P,delta,pillarq,R,rep,seed)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nick Lowery and Tristan Ursell
% 2018
%
% Structured environments fundamentally alter dynamics and stability of ecological communities
% https://www.biorxiv.org/content/early/2018/07/10/366559
% 
% This function runs a non-dimensionalized modified 2D Lotka-Volterra ecological
% simulation of two mutally killing species, and outputs images of the simulation 
% and a file containing the associated metadata. It requires the image_save.m 
% script to be in the active directory. For these simulations, the pillar spacing
% is hard coded in the initialization section rather than as an input parameter.
% Density matrices are initialized using pink noise, as fully random initialization
% leads to pathological results when competition is asymmetric.
%
% Input parameters:
% L = size of the simulation box (pixels)
% N = number of simulation iterations (dt units, specified in the script)
% D = diffusion coefficient (dimensionless)
% P = average interaction strength (dimensionless)
% delta = degree of competition asymmetry (species A gets P + delta, B gets P - delta)
% pillarq = logical (1,0); should pillars be included?
% R = pillar radius (pixels)
% rep = unique replicate identifier
% seed = sets random seet for density matrix initialization
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%*******************************************
%********* INITIALIZATION ******************
%*******************************************

% set seed
rng(seed) 

%output images?
imageq = 1;

%time step in ND units
dt = 0.01;

%local concentration discrete cutoff (will be set to 0 below this value)
dsc = 0.01;

%frequency to show results
showt = 40;

%pink noise seeding of density matrices
sigmaIC = 2;
kmax = 20;

A0 = zeros(L,L);
B0 = zeros(L,L);
X = ones(L,1)*(1:L);
Y = X';
for j = -kmax:kmax
    for k = -kmax:kmax
        if and(k ~= 0, j ~= 0)
            rand_amp = normrnd(0, sigmaIC) / sqrt(j^2 + k^2);
            theta_rand = 2 * pi * rand;
            
            A0 = A0 + sin(theta_rand) * rand_amp * sin(2 * pi / L * (j*X + k*Y)) + ...
	    	cos(theta_rand) * rand_amp * cos(2 * pi / L * (j*X + k*Y));
            
            rand_amp = normrnd(0, sigmaIC) / sqrt(j^2 + k^2);
            theta_rand = 2 * pi * rand;
            
            B0 = B0 + sin(theta_rand) * rand_amp * sin(2 * pi / L * (j*X + k*Y)) + ...
	    	cos(theta_rand) * rand_amp * cos(2 * pi / L * (j*X + k*Y));
        end
    end
    %disp(num2str(j))
end
A1 = A0;
A1(A0 < (-3 * std(A0(:)))) = -3 * std(A0(:));
A1(A0 > (3 * std(A0(:)))) = 3 * std(A0(:));
A = mat2gray(A1);

B1 = B0;
B1(B0 < (-3 * std(A0(:)))) = -3 * std(A0(:));
B1(B0 > (3 * std(A0(:)))) = 3 * std(A0(:));
B = mat2gray(B1);

A = 0.1 * A;
B = 0.1 * B;

%pillar spacing (units of R)
    % note this is a fixed lattice constant; change the coefficient to change relative pillar spacing
dx=3*R;

%diffusion convolution filter
sigma = sqrt(4*D*dt);
if mod(ceil(3*sigma), 2) == 0
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

file1 = ['2sp_LV_comp_asym_L-' num2str(L) ...
       '_N-' num2str(N) ...
       '_D-' num2str(D) ...
       '_P-' num2str(P) ...
       '_delta-' num2str(delta) ...
       '_pillars-' num2str(pillarq) ...
       '_R-' num2str(R) ...
       '_dx-' num2str(dx) ...
       '_dsc-' num2str(dsc) ...
       '_rep-' num2str(rep) ...
       '_seed-' num2str(seed)];
 
path1 = fullfile(path0, file1);
if imageq == 1
    mkdir(path1)
end

%tracking vectors (global mean, scale factor for each image)
meanA = zeros(N, 1);
meanB = zeros(N, 1);
scale_fac = zeros(N, 2);

%*******************************************
%********* PILLARS *************************
%*******************************************
if pillarq == 1
    %create spatial matrices / filters
    Xmat = ones(L,1)*(1:L);
    Ymat = Xmat';
       
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
        
    %prepare reflecting boundary conditions in grid
    conv_filt = conv2(filt_all, Gauss, 'same');
    rescale_filt = filt_all ./ conv_filt;
    rescale_filt(isnan(rescale_filt)) = 0;
else
    R = 0; dx = 0;
    filt_all = ones(L,L);
end

%calculate weight of non-pillar area for calcualting mean abundances
filt_all_weight = sum(filt_all(:));

%*******************************************
%********* SIMULATION **********************
%*******************************************

q = 0;
tic
tnow = toc;

for i = 1:N
    %diffusion convolutions
    A = conv2(A, Gauss, 'same') .* rescale_filt;
    B = conv2(B, Gauss, 'same') .* rescale_filt;
    
    %interspecies interactions
    dA = dt * A .* (1 - (A+B)) .* (1 - B / (P + delta));
    dB = dt * B .* (1 - (A+B)) .* (1 - A / (P - delta));
   
    if (sum(abs(dA(:))) + sum(abs(dB(:)))) < dsc
	t_stop = i;
	break
    end
    
    %deterministic update populations
    A = A + dA;
    B = B + dB;
 
    %hard upper bound
    A(A > 1) = 1;
    B(B > 1) = 1;
    
    %discrete cutoff (hard lower bound)
    A(A < dsc) = 0;
    B(B < dsc) = 0;
        
    %calculate tracking vectors
    if pillarq == 1
        meanA(i) = sum(A(filt_all == 1)) / filt_all_weight;
        meanB(i) = sum(B(filt_all == 1)) / filt_all_weight;
    else
        meanA(i) = mean(A(:));
        meanB(i) = mean(B(:));
    end

    scale_fac(i,1) = max(A(:));
    scale_fac(i,2) = max(B(:));
   
    %output
    if mod(i, showt) == 1
        q = q + 1;
        
        tnow_temp = toc;
        tnow = tnow_temp - tnow;
        disp(['step ' num2str(i) ' of ' num2str(N) ' for ' file1 ' in ' num2str(tnow) 's'])
        tnow = toc;
        
        %grey filter        
	grey_filt = (A == 0) .* (B == 0);

        Atemp = A;
        Btemp = B;

        if pillarq == 1
            Atemp(grey_filt == 1) = scale_fac(i,1) / 2;
            Btemp(grey_filt == 1) = scale_fac(i,2) / 2;
        end
        
        im_out(:,:,1) = uint8(Atemp / scale_fac(i,1) * 255);
        im_out(:,:,2) = uint8(Btemp / scale_fac(i,2) * 255);
        im_out(:,:,3) = uint8(Atemp / scale_fac(i,1) * 255);
        
        if imageq == 1
            image_save(im_out, fullfile(path1, [file1 '_' sprintf('%04d',q) '.png']))
        end
    end
end

%*******************************************
%********** METADATA ***********************
%*******************************************

if exist('t_stop')
    meanA_out = meanA(1:showt:t_stop - 1);
    meanB_out = meanB(1:showt:t_stop - 1);
    scale_fac_out = scale_fac(1:showt:t_stop - 1, :);
else
    meanA_out = meanA(1:showt:N);
    meanB_out = meanB(1:showt:N);
    scale_fac_out = scale_fac(1:showt:N, :);
    t_stop = NaN;
end

save([file1 '_metadata.mat'], 'L', 'N', 'D', 'P', 'delta', 'pillarq', 'R', 'dx', 'dt', 'dsc', ...
	'meanA_out', 'meanB_out', 'scale_fac_out', 'rep', 'seed', 'filt_all_weight', 't_stop')


disp(['Done with ' file1 '.'])

%exit

end
