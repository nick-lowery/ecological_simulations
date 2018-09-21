function autocorrelation_analysis(imdir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nick Lowery and Tristan Ursell
% 2018
%
% Structured environments fundamentally alter dynamics and stability of ecological communities
% https://www.biorxiv.org/content/early/2018/07/10/366559
% 
% This function calculates the spatiotemporal autocorrelation matrix for 
% simulation images generated with LV_competition_function.m, and writes it 
% to a .csv file.
%
% Input parameters:
% imdir = name of directory containing simulation images
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

drawnow

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% LOAD IMAGE STACK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filePattern = fullfile(imdir, '*.png');
files = dir(filePattern);

if isempty(files)
    error('No files match the search criteria.')
end

% extract relevant simulation parameters from file name
L = str2num(char(regexp(imdir, 'L-(\d+)', 'tokens', 'once')));
pillarq = str2num(char(regexp(imdir, 'pillars-(\d)', 'tokens', 'once')));
R = str2num(char(regexp(imdir, 'R-(\d+)', 'tokens', 'once')));
dx = str2num(char(regexp(imdir, 'dx-(\d+)', 'tokens', 'once')));
dsc = str2num(char(regexp(imdir, 'dsc-(.+)_', 'tokens', 'once')));
Nspecies = str2num(imdir(1));
nfiles = length(files);

% reorder files to correct sequence (correct for bug in original function)
[files.index] = deal(0);
for i = 1:nfiles
    filenum = str2num(char(regexp(files(i).name, '(\d+).png', 'tokens', 'once')));
    files(i).index = sprintf('%04d', filenum );
end
file_fields = fieldnames(files);
files_cell = struct2cell(files);
file_dims = size(files_cell);
files_sort = sortrows(reshape(files_cell, file_dims(1), [])',7);
ordered_files = cell2struct(reshape(files_sort', file_dims), file_fields, 1);

% load test image for size
test1 = imread(fullfile(ordered_files(1).folder, ordered_files(1).name));
[sy,sx,sz] = size(test1);

% set image subsampling
    % standard parameters write 5,000 images per simulation,
    % so subsamp = 5 returns a 1000 x 1000 autocorrelation matrix
subsamp = 5; 

% construct image arrays
all_images = zeros(sy, sx, sz, nfiles, 'uint8');
for i = 1:nfiles
    temp1 = imread(fullfile(ordered_files(i).folder, ordered_files(i).name));
    all_images(:,:,:,i) = temp1;
end
images = all_images(:,:,:,1:subsamp:end);

disp('Files imported sucessfully')

% constrain analysis to dynamic regime of simulation
    % i.e. exclude images of single species post-extinction cascade
thresh1 = 0.5 * dsc;
thresh_vec = ones(length(images),1);
for i = 2:length(images)
    diff = images(:,:,:,i) - images(:,:,:,i-1);
    thresh_vec(i) = max(abs(diff(:)));
end
no_change = find(thresh_vec < thresh1, 1, 'first');
if isempty(no_change)
    last_frame = length(images);
elseif no_change > 251
    last_frame = no_change;
else
    last_frame = 251;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% (RE)CONSTRUCT PILLAR MASK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if pillarq == 1

    % construct center-distance matrix
    Xmat = ones(sy,1)*(1:sx);
    Ymat = (1:sy)'*ones(1,sx);
    dist_mat = sqrt((Xmat-sx/2).^2+(Ymat-sy/2).^2);
    
    %form hexagonal grid position arrays
    x_cent0 = 0:dx:(L+dx);
    y_cent0 = [0:(sqrt(3)/2*dx):(L+dx*sqrt(3)/2)]+R;
    
    x_cent = zeros(length(y_cent0),length(x_cent0));
    for i = 1:length(y_cent0)
        if mod(i,2) == 1
            x_cent(i,:) = x_cent0+dx/4;
        else
            x_cent(i,:) = x_cent0-dx/4;
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
            filt_temp = ((Xmat-x_cent(j,i)).^2+(Ymat-y_cent(j,i)).^2)>=R^2;
            filt_all = filt_all.*filt_temp;
        end
    end
    
    dbound = -L/2;
    filt_all_temp = zeros(size(filt_all));
    filt_all_temp(1:round(L/2)-dbound,:) = filt_all(1:round(L/2)-dbound,:);
    filt_all_temp(round(L/2)-dbound+1:end,:) = 1;
    filt_all = filt_all_temp;
    
else
    filt_all = ones(sy,sx);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% CALCULATE SPATIOTEMPORAL AUTOCORRELATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

autocorr_mat = zeros(last_frame, last_frame);
for i = 1:last_frame
    % combine species channels into single vector, omitting pillar pixels
    tempA = double(images(:,:,1,i));
    tempB = double(images(:,:,2,i));
    tempC = double(images(:,:,3,i));
    vec1_temp = [tempA(filt_all == 1) - mean(tempA(filt_all == 1)); ...
                 tempB(filt_all == 1) - mean(tempB(filt_all == 1)); ...
                 tempC(filt_all == 1) - mean(tempC(filt_all == 1))];
    % normalize
    alpha1 = sum(vec1_temp .* vec1_temp);
    vec1 = vec1_temp ./ sqrt(alpha1);

    for j = i + 1:last_frame
        % combine species channels into single vector, omitting pillar pixels
        tempA2 = double(images(:,:,1,j));
        tempB2 = double(images(:,:,2,j));
        tempC2 = double(images(:,:,3,j));
        vec2_temp = [tempA2(filt_all == 1) - mean(tempA2(filt_all == 1)); ...
                     tempB2(filt_all == 1) - mean(tempB2(filt_all == 1)); ...
                     tempC2(filt_all == 1) - mean(tempC2(filt_all == 1))];
        % normalize
        alpha2 = sum(vec2_temp .* vec2_temp);
        vec2 = vec2_temp ./ sqrt(alpha2);
        
        % calculate correlation coefficient
        autocorr_mat(i,j) = sum(vec1 .* vec2);
    end

    disp(['Finished with ' num2str(i) ' of ' num2str(last_frame)])
end

% mirror upper triangle to fill in matrix, add correlation of 1 on diagonal
autocorr_mat = autocorr_mat + autocorr_mat' + eye(last_frame);

% write file to disk
csvwrite([imdir '_acf.csv'], autocorr_mat)

disp(['Files written - Finished!.'])

end


