function LV_competition_two_species_asym_pink_lattice(L,N,D,P,delta,pillarq,R,rep,seed)

%Tristan & Nick -- LV Spatial Competition Model
%Sept 2017
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
%N=1000;

%time step in ND units
dt=0.01;

%ND diffusion coefficient
%D=15;

%interaction parameter (in units of carrying capacity)
%P=Inf;
%P=0.1;
%P=[-10,-1,-0.1,-0.01,0.01,0.1,1,10];

%local concentration discrete cutoff (will be set to 0 below this value)
dsc=0.001;

%intialize the density matrices
    %'drop' in center
%    ini_mat=zeros(L,L);
%    ini_mat(round(L/2),round(L/2))=1;
%    ini_mat=double(bwdist(ini_mat)<250);
%    A=1/3*rand(L,L);
%    B=1/3*rand(L,L);
%    A(A>B)=1;
%    A(B>=A)=0;
%    B(A>B)=0;
%    B(B>=A)=1;
%    A=A.*ini_mat;
%    B=B.*ini_mat;

%uniform random seeding
%A=1/10*rand(L,L);
%B=1/10*rand(L,L);

%{
%sparse seeding
initA = rand(L,L);
initB = rand(L,L);
A = rand(L,L);
B = rand(L,L);
A(initA > 0.01) = 0;
B(initB > 0.01) = 0;
%}

%%{
%pink noise seeding
sigmaIC=2;
kmax=20;

A0=zeros(L,L);
B0=zeros(L,L);
X=ones(L,1)*(1:L);
Y=X';
for j=-kmax:kmax
    for k=-kmax:kmax
        if and(k~=0,j~=0)
            rand_amp=normrnd(0,sigmaIC)/sqrt(j^2+k^2);
            theta_rand=2*pi*rand;
            
            A0=A0+sin(theta_rand)*rand_amp*sin(2*pi/L*(j*X+k*Y))+cos(theta_rand)*rand_amp*cos(2*pi/L*(j*X+k*Y));
            
            rand_amp=normrnd(0,sigmaIC)/sqrt(j^2+k^2);
            theta_rand=2*pi*rand;
            
            B0=B0+sin(theta_rand)*rand_amp*sin(2*pi/L*(j*X+k*Y))+cos(theta_rand)*rand_amp*cos(2*pi/L*(j*X+k*Y));
        end
    end
    disp(num2str(j))
end
A1=A0;
A1(A0<(-3*std(A0(:))))=-3*std(A0(:));
A1(A0>(3*std(A0(:))))=3*std(A0(:));
A=mat2gray(A1);

B1=B0;
B1(B0<(-3*std(A0(:))))=-3*std(A0(:));
B1(B0>(3*std(A0(:))))=3*std(A0(:));
B=mat2gray(B1);

A=0.1*A;
B=0.1*B;
%}

%frequency to show results
showt=40;

%use pillars?
%pillarq=1;

%pillar radius
%R=15;
%R=[5, 10, 15, 20, 25, 30, 35, 50];

%pillar spacing (units of R)
%dx_mult=4;
dx=3*R;
%dx=[2.5, 3, 3.5, 4, 4.5, 5]*R;


%*************************************************************************

%diffusion convolution filter
sigma=sqrt(4*D*dt);
if mod(ceil(4*sigma),2)==0
    win_sz=ceil(4*sigma)+1;
else
    win_sz=ceil(4*sigma);
end
Gauss=fspecial('gaussian',[win_sz,win_sz],sigma);
rescale_filt=1./conv2(ones(size(A)),Gauss,'same');

if L<win_sz
    warning('Diffusion window size is larger than simulation box.')
end

%file IO
path0=cd;

file1=['2sp_LV_comp_asym_L-' num2str(L) ...
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
 
path1=fullfile(path0,file1);
if imageq==1
    mkdir(path1)
end

%tracking vectors (global mean, point, scale factor for each image)
meanA=zeros(N,1);
meanB=zeros(N,1);
scale_fac=zeros(N,2);

%*******************************************
%********* PILLARS *************************
%*******************************************
if pillarq==1
    %create spatial matrices / filters
    Xmat=ones(L,1)*(1:L);
    Ymat=Xmat';
      
    %form hexagonal grid position arrays
    x_cent0=0:dx:(L+dx);
    y_cent0=[0:(sqrt(3)/2*dx):(L+dx*sqrt(3)/2)]+R;
    
    x_cent=zeros(length(y_cent0),length(x_cent0));
    for i=1:length(y_cent0)
        if mod(i,2)==1
            x_cent(i,:)=x_cent0+dx/4;
        else
            x_cent(i,:)=x_cent0-dx/4;
        end
    end
    
    y_cent=zeros(size(x_cent));
    for j=1:length(x_cent0)
        y_cent(:,j)=y_cent0';
    end
    
    %create pillars
    filt_all=ones(size(Xmat));
    for i=1:length(x_cent0)
        for j=1:length(y_cent0)
            filt_temp=((Xmat-x_cent(j,i)).^2+(Ymat-y_cent(j,i)).^2)>=R^2;
            filt_all=filt_all.*filt_temp;
        end
    end
    
    dbound=-L/2;
    filt_all_temp=zeros(size(filt_all));
    filt_all_temp(1:round(L/2)-dbound,:)=filt_all(1:round(L/2)-dbound,:);
    %filt_all_temp(round(L/2)-dbound+1:end,:)=filt_all2(round(L/2)-dbound+1:end,:);
    filt_all_temp(round(L/2)-dbound+1:end,:)=1;
    filt_all=filt_all_temp;
    
    
    %sum weights of active areas for calcualting mean abundances
    filt_all_weight=sum(filt_all(:));
    
    %remove objects on edge
    %filt_all=double(~imclearborder(~filt_all));
    
    %prepare reflecting boundary conditions in grid
    conv_filt=conv2(filt_all,Gauss,'same');
    rescale_filt=filt_all./conv_filt;
    rescale_filt(isnan(rescale_filt))=0;
else
    filt_all = ones(L,L);
    filt_all_weight = sum(filt_all(:));
end
%*******************************************
%*********  END PILLARS ********************
%*******************************************


%{
%*******************************************
%Setup calculation of structure factor
%*******************************************

%generate center-distance matrix
temp1=zeros(L,L);
temp1(round(L/2),round(L/2))=1;
dist_mat1=bwdist(temp1);

%generate center-distance vector
dist_vec=dist_mat1(:);

bins1=0:2:L;
for j=1:length(bins1)-1
    temp2=find(and(dist_vec>bins1(j),dist_vec<=bins1(j+1)));
    dist_els(j).vec=temp2;
    dist_els(j).length=length(temp2);
end

struc_vec=zeros(floor(N/showt),length(bins1)-1);
%*******************************************
%}

q=0;
if plotq
    figure;
end
tic
tnow=toc;
for i=1:N
    %convolutions
    A=conv2(A,Gauss,'same').*rescale_filt;
    B=conv2(B,Gauss,'same').*rescale_filt;
    
    %population networks    
    descript1='simple_intrans_comp';
    dA=dt*A.*(1-(A+B)).*(1-B/(P+delta));
    dB=dt*B.*(1-(A+B)).*(1-A/(P-delta));
   
    if (sum(abs(dA(:))) + sum(abs(dB(:)))) < dsc
	t_stop = i;
	break
    end
    
    %deterministic update populations
    A=A+dA;
    B=B+dB;
 
    %hard upper bound
    A(A>1)=1;
    B(B>1)=1;
    
    %discrete cutoff (hard lower bound)
    A(A<dsc)=0;
    B(B<dsc)=0;
        
    %calculate tracking vectors
    if pillarq==1
        meanA(i)=sum(A(filt_all==1))/filt_all_weight;
        meanB(i)=sum(B(filt_all==1))/filt_all_weight;
    else
        meanA(i)=mean(A(:));
        meanB(i)=mean(B(:));
    end

    scale_fac(i,1)=max(A(:));
    scale_fac(i,2)=max(B(:));
   
    %output
    if mod(i,showt)==1
        q=q+1;
        
        tnow_temp=toc;
        tnow=tnow_temp-tnow;
        disp(['step ' num2str(i) ' of ' num2str(N) ' for ' file1 ' in ' num2str(tnow) 's'])
        tnow=toc;
        
        %grey filter        
	grey_filt=(A==0).*(B==0);

        Atemp=A;
        Btemp=B;

        if pillarq==1
            Atemp(grey_filt==1)=scale_fac(i,1)/2;
            Btemp(grey_filt==1)=scale_fac(i,2)/2;
        end
        
        im_out(:,:,1)=uint8(Atemp/scale_fac(i,1)*255);
        im_out(:,:,2)=uint8(Btemp/scale_fac(i,2)*255);
        im_out(:,:,3)=uint8(Atemp/scale_fac(i,1)*255);
        
        if plotq==1
            %imagesctsu(im_out,[0 1])
            imagesctsu(im_out)
            xlabel('X')
            ylabel('Y')
            title(['time = ' num2str(i)])
            drawnow
        end
        
        if imageq==1            
            %im_out2=uint8(im_out*255);
            
            image_save(im_out,fullfile(path1,[file1 '_' sprintf('%04d',q) '.png']))
        end
        

        %{
        %calculate structure factor
        struc_mat=abs(fftshift(fft2(Atemp)));
        
        %rotationally average
        struc_temp=struc_mat(:);
        
        for j=1:length(bins1)-1
            struc_vec(q,j)=mean(struc_temp(dist_els(j).vec));
        end
        %}
    end
end

%***********************************************
%************ FIGURES **************************
%***********************************************
%create time vector (in units of doubling time)
%{
dt_dt=2;
t_vec=(1:dt_dt:i-1)*dt;

h1=figure('rend','painters','pos',[200 200 1100 500]);
hold on
plot(t_vec,meanA(1:dt_dt:i-1),'r','linewidth',2)
plot(t_vec,meanB(1:dt_dt:i-1),'g','linewidth',2)
ylabel('System-wide Abundance (carrying capacity)')
xlabel('Time (doubling periods)')
title(['D = ' num2str(D) ', R = ' num2str(R) ', dx = ' num2str(dx)])
xlim([0 max(t_vec)])
box on

saveas(h1,[file1 '_trace.fig'])
%}

if exist('t_stop')
    meanA_out = meanA(1:showt:t_stop-1);
    meanB_out = meanB(1:showt:t_stop-1);
    scale_fac_out = scale_fac(1:showt:t_stop-1, :);
else
    meanA_out = meanA(1:showt:N);
    meanB_out = meanB(1:showt:N);
    scale_fac_out = scale_fac(1:showt:N, :);
    t_stop = NaN;
end

save([file1 '_metadata.mat'], 'L', 'N', 'D', 'P', 'delta', 'pillarq', 'R', 'dx', 'dt', 'dsc', ...
	'meanA_out', 'meanB_out', 'scale_fac_out', 'rep', 'seed', 'filt_all_weight', 't_stop')


disp(['Done with ' file1 '.'])

exit

end










