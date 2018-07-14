%%% Reconstruct pillar mask %%%%%%%%%%%%%%

% define relevant parameters
D = 15;
R = 30;
dx = 5*R;
sy = 500;
sx = 500;
L = 500;

var_R = 0;  % max ~ 0.5
var_dx = 0; % relative to R; depends on dx

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

%%% jitter pillar placement

% select random angle
theta = 360*rand(size(x_cent));

% generate random displacement within +/- var_dx
d = var_dx .* rand(size(x_cent)) * R; 

x_jitt = x_cent + d .* cosd(theta);
y_jitt = y_cent + d .* sind(theta);

%create pillars
filt_all = ones(size(Xmat));
for i = 1:length(x_cent0)
    for j = 1:length(y_cent0)
        R_jitter = (-var_R + 2*var_R .* rand(1))*R + R;
        if R_jitter <= D/2
            R_jitter = D/2 + 0.1;
        end
        filt_temp = ((Xmat-x_jitt(j,i)).^2+(Ymat-y_jitt(j,i)).^2) >= R_jitter^2;
        filt_all = filt_all.*filt_temp;
    end
end

dbound = -L/2;
filt_all_temp = zeros(size(filt_all));
filt_all_temp(1:round(L/2)-dbound,:) = filt_all(1:round(L/2)-dbound,:);
filt_all_temp(round(L/2)-dbound+1:end,:) = 1;
filt_all = filt_all_temp;

%%% End pillar mask %%%%%%%%%%%%%%%%%%

imagesc(filt_all)

%{
grids = dir('jitter*.mat');
grid_names = {grids(:).name};
str = ['R-' num2str(R) '_dx-' num2str(dx/R) '_vR-' num2str(jitt_R) '_vdx-' num2str(jitt_dx) '_num-' num2str(jitt_num)];
filt_all_temp = load(string(grid_names(find(contains(grid_names, str)))));
filt_all = filt_all_temp.filt_all;
%}