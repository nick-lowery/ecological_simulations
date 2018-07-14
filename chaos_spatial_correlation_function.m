function chaos_spatial_correlation_function(imdir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nick Lowery and Tristan Ursell
% 2018
%
% Structured environments fundamentally alter dynamics and stability of ecological communities
% https://www.biorxiv.org/content/early/2018/07/10/366559
% 
% This function calculates the spatiotemporal correlations between simulation
% images generated with the LV_competition_chaos function, and writes a data
% file containing the temportal correlation vectors.
%
% Input parameters:
% imdir = name of directories containing simulations to be compared,
%	as written in the LV_competition_chaos.m script, and truncated to 
% 	omit replicates
% For example, for simulation directories
% 3sp_LV_comp_chaos...dinit-0.1_rep-1
% 3sp_LV_comp_chaos...dinit-0.1_rep-2
% 3sp_LV_comp_chaos...dinit-0.1_rep-3
% ...
% imdir = 3sp_LV_comp_chaos...dinit-0.1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

drawnow

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% FIND DIRECTORIES AND PARAMETERS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find directories that match imdir
nreps = 10;
dirnames = squeeze(strings(nreps,1));
for i = 1:nreps
    dirnames(i) = [imdir '_rep-' num2str(i)];
end

% extract relevant simulation parameters
L = str2num(char(regexp(imdir, 'L-(\d+)', 'tokens', 'once')));
pillarq = str2num(char(regexp(imdir, 'pillars-(\d)', 'tokens', 'once')));
R = str2num(char(regexp(imdir, 'R-(\d+)', 'tokens', 'once')));
dx = str2num(char(regexp(imdir, 'dx-(\d+)', 'tokens', 'once')));
dinit = str2num(char(regexp(imdir, 'dinit-(\d.\d+)', 'tokens', 'once')));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% PILLARS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sx = L;
sy = L;

if pillarq == 1
    %%% Reconstruct pillar mask %%%%%%%%%%%%%%

    % construct center-distance matrix
    Xmat = ones(sy,1)*(1:sx);
    Ymat = (1:sy)'*ones(1,sx);
    dist_mat = sqrt((Xmat-sx/2).^2 + (Ymat-sy/2).^2);

    %form hexagonal grid position arrays
    x_cent0 = 0:dx:(L+dx);
    y_cent0 = [0:(sqrt(3)/2*dx):(L+dx*sqrt(3)/2)] + R;

    x_cent = zeros(length(y_cent0),length(x_cent0));
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
            filt_temp = ((Xmat-x_cent(j,i)).^2 + (Ymat-y_cent(j,i)).^2) >= R^2;
            filt_all = filt_all .* filt_temp;
        end
    end

    dbound = -L/2;
    filt_all_temp = zeros(size(filt_all));
    filt_all_temp(1:round(L/2) - dbound,:) = filt_all(1:round(L/2) - dbound,:);
    filt_all_temp(round(L/2) - dbound + 1:end,:) = 1;
    filt_all = filt_all_temp;

    %%% End pillar mask %%%%%%%%%%%%%%%%%%
else
    filt_all = ones(L);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% CORRELATIONS %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set image subsampling frequency (for speed)
	% note total number of images will depend on N and showt of simulation
subsamp = 5;
img_vec = 1:subsamp:5000;
nfiles = length(img_vec);

% initialize correlation and image stack matrices
ncomps = round(nreps * (nreps - 1)/2);
cor_mat = zeros(nfiles, ncomps);
images1 = zeros(sx, sy, 3, nfiles, 'uint8');
images2 = zeros(sx, sy, 3, nfiles, 'uint8');
compcount = 1;

disp('Starting correlation calculations...')

% loop through stacks for pairwise correlations
for i = 1:length(dirnames)
    for j = i + 1:length(dirnames)
        
        % read filenames
        filestack1 = dir(fullfile([char(dirnames(i)) '/*.png']));
        filestack2 = dir(fullfile([char(dirnames(j)) '/*.png']));
        
        % read in image stacks
        for k = 1:nfiles
            temp1 = imread(fullfile(filestack1(img_vec(k)).folder, ...
                                    filestack1(img_vec(k)).name));
            temp2 = imread(fullfile(filestack2(img_vec(k)).folder, ...
                                    filestack2(img_vec(k)).name));
            images1(:,:,:,k) = temp1;
            images2(:,:,:,k) = temp2;
        end

	% load scale factors
	scale_temp_1 = load(strcat(dirnames(i), '_metadata.mat'), 'scale_fac_out');
	scale_fac_1 = scale_temp_1.scale_fac_out;
	scale_fac_1 = scale_fac_1(img_vec, :);
        
	scale_temp_2 = load(strcat(dirnames(j), '_metadata.mat'), 'scale_fac_out');
	scale_fac_2 = scale_temp_2.scale_fac_out;
	scale_fac_2 = scale_fac_2(img_vec, :);

        % calculate correlation
        for m = 1:length(images1)
	    % separate image channels
            tempA1 = double(images1(:,:,1,m));
            tempB1 = double(images1(:,:,2,m));
            tempC1 = double(images1(:,:,3,m));
	    
	    % rescale channels
	    tempA1 = tempA1 .* scale_fac_1(m,1);
	    tempB1 = tempB1 .* scale_fac_1(m,2);
	    tempC1 = tempC1 .* scale_fac_1(m,3);

	    % vectorize and normalize image
            vec1_temp = [tempA1(filt_all == 1) - mean(tempA1(filt_all == 1)); ...
                         tempB1(filt_all == 1) - mean(tempB1(filt_all == 1)); ...
                         tempC1(filt_all == 1) - mean(tempC1(filt_all == 1))];
            alpha1 = sum(vec1_temp .* vec1_temp);
            vec1 = vec1_temp ./ sqrt(alpha1);
            
	    % separate image channels
            tempA2 = double(images2(:,:,1,m));
            tempB2 = double(images2(:,:,2,m));
            tempC2 = double(images2(:,:,3,m));
	    
	    % rescale channels
	    tempA2 = tempA2 .* scale_fac_2(m,1);
	    tempB2 = tempB2 .* scale_fac_2(m,2);
	    tempC2 = tempC2 .* scale_fac_2(m,3);
            
	    % vectorize and normalize image
	    vec2_temp = [tempA2(filt_all == 1) - mean(tempA2(filt_all == 1)); ...
                         tempB2(filt_all == 1) - mean(tempB2(filt_all == 1)); ...
                         tempC2(filt_all == 1) - mean(tempC2(filt_all == 1))];
            alpha2 = sum(vec2_temp .* vec2_temp);
            vec2 = vec2_temp ./ sqrt(alpha2);
            
	    % calculate correlation coefficient
            cor_mat(m,compcount) = sum(vec1 .* vec2);
        end        
        disp(['Finished correlation of rep ' num2str(i) ' with rep ' num2str(j)])
        compcount = compcount + 1;
    end
end

cor_mean = mean(cor_mat, 2);
cor_sd = std(cor_mat, 0, 2);

save([imdir '_cor.mat'], 'cor_mat', 'cor_mean', 'cor_sd', 'L', 'R', 'dx', 'dinit')

disp(['Files written, done with ' imdir])

%exit

end
