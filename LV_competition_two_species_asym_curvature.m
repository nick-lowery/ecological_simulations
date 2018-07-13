function LV_competition_two_species_asym_curvature_noclock(W,H,D,P,delta,R)

%Tristan & Nick -- LV Spatial Competition Model
%Sept 2017
% ND = non-dimensional
%

%output images
imageq=1;

%box size (grid points)
%H=220;
%W=350;

%number of time steps
%N=100000;

%time step in ND units
dt=0.01;

%ND diffusion coefficient
%D=15;

%interaction parameter (in units of carrying capacity)
%P=0.1;

%competition asymmetry
%delta=0.024;

%local concentration discrete cutoff (will be set to 0 below this value)
dsc=0.001;

%intialize the density matrices
%%{
    %block seeding
    A = zeros(H,W);
    A(:,1:W/2) = 0.3;
    B = zeros(H,W);
    B(:,W/2:W) = 0.3;
%}

%frequency to show results
showt=50;

%pillar radius
%R=50;

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

if min(H,W)<win_sz
    warning('Diffusion window size is larger than simulation box.')
end

%file IO
file1=['2sp_LV_comp_curve_H-' num2str(H) ...
       '_W-' num2str(W) ...
       '_D-' num2str(D) ...
       '_P-' num2str(P) ...
       '_delta-' num2str(delta) ...
       '_R-' num2str(R) ...
       '_dsc-' num2str(dsc)];
 
%*******************************************
%********* PILLARS *************************
%*******************************************

%create spatial matrices / filters
Xmat=ones(H,1)*(1:W);
Ymat=(ones(W,1)*(1:H))';

%pillar locations
x_cent=[W/2 W/2];
y_cent=[R H-R];    

%create pillars
filt_all=ones(size(Xmat));
for i=1:length(x_cent)
    filt_temp=((Xmat-x_cent(i)).^2+(Ymat-y_cent(i)).^2)>=R^2;
    filt_all=filt_all.*filt_temp;
end

%prepare reflecting boundary conditions in grid
conv_filt=conv2(filt_all,Gauss,'same');
rescale_filt=filt_all./conv_filt;
rescale_filt(isnan(rescale_filt))=0;

%*******************************************
%*********  END PILLARS ********************
%*******************************************

q=0;
arrest = 0;
tic
tnow=toc;
while arrest < 1
    %track image for comparison
    Aprev=A;
    Bprev=B;
    
    %convolutions
    A=conv2(A,Gauss,'same').*rescale_filt;
    B=conv2(B,Gauss,'same').*rescale_filt;
    
    %population networks    
    dA=dt*A.*(1-(A+B)).*(1-B/(P+delta));
    dB=dt*B.*(1-(A+B)).*(1-A/(P-delta));
    
    %deterministic update populations
    A=A+dA;
    B=B+dB;
    
    %hard upper bound
    A(A>1)=1;
    B(B>1)=1;
    
    %discrete cutoff (hard lower bound)
    A(A<dsc)=0;
    B(B<dsc)=0;
    
    if and(max(abs(Aprev-A))<(100*eps), max(abs(Bprev-B))<(100*eps))
        disp('Simulation terminated due to extinction or inactivity.')
        arrest = 1;
	break
    end
    
    q = q+1;

    %output
    if mod(q,showt)==1
        tnow_temp=toc;
        tnow=tnow_temp-tnow;
        disp(['step ' num2str(i) ' for ' file1 ' in ' num2str(tnow) 's'])
        tnow=toc;
    end
end

%***********************************************
%************ OUTPUT ***************************
%***********************************************

%grey filter 
Atemp=A;
Btemp=B;
Atemp(filt_all==0)=1/2;
Btemp(filt_all==0)=1/2;

im_out = zeros(H,W,3,'uint8');
im_out(:,:,1)=uint8(Atemp*255);
im_out(:,:,2)=uint8(Btemp*255);
im_out(:,:,3)=uint8(Atemp*255);

%***********************************************
%************ CURVATURE ************************
%***********************************************

% hard boundary, based on A majority
im_bin = A > B;

if sum(im_bin(:)) == sum(filt_all(:))
    K = NaN;
    extinct = 1;
else
    extinct = 0;

    % extract curve
    im_curve = edge(im_bin, 'Sobel') - edge(filt_all, 'Sobel');
    curve = im_curve == 1;

    % define curve coordinates
    curve_props = regionprops(curve, 'Extrema');

    % top right corner
    x1 = curve_props.Extrema(2,1);
    y1 = curve_props.Extrema(2,2);
    % right side of midpoint
    x2 = max(curve_props.Extrema(:,1));
    y2 = 0.5*(max(curve_props.Extrema(:,2))+min(curve_props.Extrema(:,2)));
    % bottom right corner
    x3 = curve_props.Extrema(5,1);
    y3 = curve_props.Extrema(5,2);

    % calculate curvature
    a = sqrt((x1-x2)^2+(y1-y2)^2);
    b = sqrt((x2-x3)^2+(y2-y3)^2);
    c = sqrt((x3-x1)^2+(y3-y1)^2);
    A = 1/2*abs((x1-x2)*(y3-y2)-(y1-y2)*(x3-x2)); 
    if a*b*c == 0
	K = 0;
    else 
        K = 4*A/(a*b*c); 
    end
end

disp(['extinction = ' num2str(extinct) ', K = ' num2str(K)])

%***********************************************
%************ FIGURES **************************
%***********************************************

if imageq==1            
    image_save(im_out,[file1 '.png'])
    save([file1 '.mat'], 'W','H','dt','dsc','D','P','delta','R','K','extinct')
end

disp(['Done with ' file1 '.'])

exit

end








