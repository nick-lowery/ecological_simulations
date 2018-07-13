function LV_competition_function_jitter_grid(L,N,D,P,pillarq,R,dx,rep,jitt_R,jitt_dx,jitt_num,seed)

% Tristan & Nick -- LV Spatial Competition Model
% Sept 2017
% ND = non-dimensional
%

%rng(round(exp(pi/2)*1e8))
rng(seed)

%show plots
plotq=0;

%output images
imageq=1;

%box size (grid points)
%L=500;

%number of time steps
%N=100000;

%time step in ND units
dt=0.01;

%ND diffusion coefficient
%D=25;

%Interaction parameter (units of carrying capacity)
%P=1

%carrying capacity (units of critical conc)
%K=100;

%local concentration discrete cutoff (will be set to 0 below this value)
dsc=0.01;

%intialize the density matrices
A=1/10*rand(L,L);
B=1/10*rand(L,L);
C=1/10*rand(L,L);

%frequency to show results
showt=20;

%*************************************************************************

%diffusion convolution filter
sigma=sqrt(4*D*dt);
if mod(ceil(3*sigma),2)==0
    win_sz=ceil(3*sigma)+1;
else
    win_sz=ceil(3*sigma);
end
Gauss=fspecial('gaussian',[win_sz,win_sz],sigma);
rescale_filt=1./conv2(ones(size(A)),Gauss,'same');

if L<win_sz
    warning('Diffusion window size is larger than simulation box.')
end

%file IO
if imageq==1
    path0=cd;
    
    file1=['3sp_LV_comp_jitter_grid_L-' num2str(L) ...
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
    
    path1=fullfile(path0,file1);
    mkdir(path1)
    
    q=0;
end

%tracking vectors (global mean, point, scale factor for each image)
meanA=zeros(N,1);
meanB=zeros(N,1);
meanC=zeros(N,1);

scale_fac=zeros(N,3);

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
filt_all_weight=sum(filt_all(:));
    
%prepare reflecting boundary conditions in grid
conv_filt=conv2(filt_all,Gauss,'same');
rescale_filt=filt_all./conv_filt;
rescale_filt(isnan(rescale_filt))=0;

%*******************************************
%*********  END PILLARS ********************
%*******************************************

if plotq
    figure;
end
for i=1:N
    %convolutions
    A=conv2(A,Gauss,'same').*rescale_filt;
    B=conv2(B,Gauss,'same').*rescale_filt;
    C=conv2(C,Gauss,'same').*rescale_filt;
        
    %population networks    
    %%{
    descript1='3sp_simple_intrans_comp';
    dA=dt*A.*(1-(A+B+C)).*(1-C/P);
    dB=dt*B.*(1-(A+B+C)).*(1-A/P);
    dC=dt*C.*(1-(A+B+C)).*(1-B/P);
    %}
    
    %deterministic update populations
    A=A+dA;
    B=B+dB;
    C=C+dC;
    
    %hard upper bound
    A(A>1)=1;
    B(B>1)=1;
    C(C>1)=1;
    
    %discrete cutoff (hard lower bound)
    A(A<dsc)=0;
    B(B<dsc)=0;
    C(C<dsc)=0;
    
    %calculate tracking vectors
    if pillarq==1
        meanA(i)=sum(A(filt_all==1))/filt_all_weight;
        meanB(i)=sum(B(filt_all==1))/filt_all_weight;
        meanC(i)=sum(C(filt_all==1))/filt_all_weight;
    else
        meanA(i)=mean(A(:));
        meanB(i)=mean(B(:));
        meanC(i)=mean(C(:));
    end
    
    scale_fac(i,1)=max(A(:));
    scale_fac(i,2)=max(B(:));
    scale_fac(i,3)=max(C(:));
%    scale_fac(i)=max([A(:);B(:);C(:)]);
    
    %output
    if mod(i,showt)==1
        disp(['step ' num2str(i) ' of ' num2str(N) ' for ' file1])
	drawnow('update')        
 
        %grey filter
        grey_filt=(A==0).*(B==0).*(C==0);
        
        Atemp=A; 
        Btemp=B;
        Ctemp=C;
           
        if pillarq==1
%             Atemp(grey_filt==1)=scale_fac(i)/2;
%             Btemp(grey_filt==1)=scale_fac(i)/2;
%             Ctemp(grey_filt==1)=scale_fac(i)/2;
            Atemp(grey_filt==1)=scale_fac(i,1)/2;
            Btemp(grey_filt==1)=scale_fac(i,2)/2;
            Ctemp(grey_filt==1)=scale_fac(i,3)/2;
        end
        
        im_out(:,:,1)=uint8(Atemp/scale_fac(i,1)*255);
        im_out(:,:,2)=uint8(Btemp/scale_fac(i,2)*255);
        im_out(:,:,3)=uint8(Ctemp/scale_fac(i,3)*255);

%         im_out(:,:,1)=Atemp;
%         im_out(:,:,2)=Btemp;
%         im_out(:,:,3)=Ctemp;
        
        if plotq==1
            imagesctsu(im_out)
            xlabel('X')
            ylabel('Y')
            title(['time = ' num2str(i)])
            drawnow
        end
        
        if imageq==1
            q=q+1;
            
%             im_out2=uint8(im_out/scale_fac(i)*255);
%             
%             image_save(im_out2,fullfile(path1,['LV_comp ' file1 '_' descript1 '_D-' num2str(D) '_R-' num2str(R) '_dx-' num2str(dx)...
%                 '_dsc-' num2str(dsc) '_' num2str(q) '.png']))
            
            image_save(im_out,fullfile(path1,[file1 '_' sprintf('%04d',q) '.png']))
        end
    end
end

%***********************************************
%************ FIGURES **************************
%***********************************************

%{
%create time vector (in units of doubling time)
dt_dt=2;
t_vec=(1:dt_dt:i-1)*dt;

h1=figure('rend','painters','pos',[200 200 1100 500]);
subplot(1,2,1)
hold on
plot(t_vec,meanA(1:dt_dt:i-1),'r','linewidth',2)
plot(t_vec,meanB(1:dt_dt:i-1),'g','linewidth',2)
plot(t_vec,meanC(1:dt_dt:i-1),'b','linewidth',2)
ylabel('System-wide Abundance (carrying capacity)')
xlabel('Time (doubling periods)')
title(['D = ' num2str(D) ', R = ' num2str(R) ', dx = ' num2str(dx)])
xlim([0 max(t_vec)])
box on

subplot(1,2,2)
hold on
plot(t_vec,specA(1:dt_dt:i-1),'r','linewidth',2)
plot(t_vec,specB(1:dt_dt:i-1),'g','linewidth',2)
plot(t_vec,specC(1:dt_dt:i-1),'b','linewidth',2)
ylabel('Fixed Point Abundance (carrying capacity)')
xlabel('Time (doubling periods)')
xlim([0 max(t_vec)])
box on

saveas(h1,[file1 '_trace.fig'])

h2=figure('rend','painters','pos',[200 200 600 500]);
[NA,XA]=hist(meanA(35:i-1),30);
[NB,XB]=hist(meanB(35:i-1),30);
[NC,XC]=hist(meanC(35:i-1),30);
hold on
bar(XA,NA/(sum(NA)*(XA(2)-XA(1))),1,'r','linestyle','none')
bar(XB,NB/(sum(NB)*(XB(2)-XB(1))),1,'g','linestyle','none')
bar(XC,NC/(sum(NC)*(XC(2)-XC(1))),1,'b','linestyle','none')
xlabel('A (red), B (green), C (blue) [carrying capacity]')
ylabel('Abundance Probability Density')
box on

saveas(h2,[file1 '_histogram.fig'])

%h3=figure('rend','painters','pos',[200 200 1100 500]);
h3=figure('rend','painters','pos',[200 200 600 500]);
cmap1=jet(length(1:i-1));

%subplot(1,2,1)
hold on
for j=1:10:i-1
    plot3(meanA(j),meanB(j),meanC(j),'.','color',cmap1(j,:))
    disp(num2str(j))
end

max_mean=max([meanA(:);meanB(:);meanC(:)]);

%plot3(meanA(1:i-1)/K,meanB(1:i-1)/K,meanC(1:i-1)/K,'.','color',[0 0.3 1])
xlabel('A')
ylabel('B')
zlabel('C')
xlim([0 max_mean])
ylim([0 max_mean])
zlim([0 max_mean])
view([135,25])
title(['System-wide Abundance (carrying capacity)'])
box on

%{
subplot(1,2,2)
hold on
for j=1:5:i-1
    plot3(specA(j)/K,specB(j)/K,specC(j)/K,'.','color',cmap1(j,:))
    disp(num2str(j))
end

%plot3(specA(1:i-1)/K,specB(1:i-1)/K,specC(1:i-1)/K)
xlabel('A')
ylabel('B')
zlabel('C')
xlim([0 1])
ylim([0 1])
zlim([0 1])
view([135,25])
title(['Fixed Point Abundance (carrying capacity)'])
box on
%}

%saveas(h3,[file1 '_phase_plot.fig'])

%}

meanA_out = meanA(1:showt:N);
meanB_out = meanB(1:showt:N);
meanC_out = meanC(1:showt:N);

scale_fac_out = scale_fac(1:showt:N, :)

save([file1 '_metadata.mat'], 'L', 'N', 'D', 'P', 'pillarq', 'R', 'dx', 'dt', 'dsc', ...
	'meanA_out', 'meanB_out', 'meanC_out', 'scale_fac_out', 'rep', 'seed', ...
	'jitt_R', 'jitt_dx', 'jitt_num', 'filt_all_weight')

disp(['Done with ' file1 '.'])

exit

end
