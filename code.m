

openfilepath= 'D:\Matlab_files\Alfazal\test2';
% % fname='test2';


x_samples=100 ;  % sample numbers in the x scan
x_step=0.031 ;    % step length in the x scan direction, mm
y_samples=100 ;  % sample numbers in the y scan
y_step=0.031 ;    % step length in the y scan direction, mm
z_points=2000 ;  % sample numbers of A-Scan
D_trans=3 ;      % diameter of the transducer, mm
f_trans=6 ;      % the focal length of the focused      transducer, mm
 
rate_sample=100 ; % sample rate, Mhz
v_sound=1.5 ;          % velocity of sound, mm/us
% originaldata=(originaldata-540)/2048;
% openfilename=[openfilepath,fname,'.dat'];  % file's path and name
openfilename=[openfilepath,'.dat'];  % file's path and name
fid=fopen(openfilename,'r');          % read the data into fid from the file
AScans=fscanf(fid,'%f');
I=reshape(AScans,z_points,x_samples,y_samples); % transmit the 1-D original data into 3-D arry with 
                                                      % x=z_points,
                                                      % y=x_samples, 
                                                      % z=y_samples.
                                                                                   
%%%%%%%%%%%%%%%%%%%%%%%% 3D reconstruction process %%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%original data analysis%%%%%%%%%%%%%%%%%%%%%%%
I_shift=shiftdim(I,3); % shift the dim, 1500*200*71 --(3)> 1500*200*71 (x_samples=200,y_samples=71, z_points=1500)
n_section_l=300;       % the efficient start position take from the original signal
n_section_h=550;      % the effficient end position take from the original signal
I_use=(I(n_section_l:n_section_h,:,:));  % use the data between the range 
                                             % of n_section_l<=x<=n_section_h 
for ix=1:x_samples
    for iy=1:y_samples
        blf=hilbert(I_use(:,ix,iy));              % discrete Fourie transform
        def=abs(blf);                  % mode of blf
        I_use(:,ix,iy)=def;
    end
end
% for iy=1:y_samples
%     I_use_sort=sort(I_use(:,:,i00y),'descend');
%     I_use_maximum=max(max(sum(I_use_sort(1:20,:,iy))))
% end

I_max = max(max(max(I_use))) ;
I_ave=mean(I_use);
I_ave_max=max(max(I_max)) ;

% Ele_error=find(I_use>I_ave_max*10);
% I_use(Ele_error)=I_ave_max/10;

l_data=n_section_h-n_section_l;
I_use=I_use/I_ave_max;     % normalization

%%%%%%%%%%%%%%%%%%%%%%
% x=(1:x_samples)*x_step;
% z=(1:l_data)/rate_sample*v_sound;
% section=I_use(:,:,5); % get the section data of z_real*x_real at y_real=i
% figure
% imshow(section,'XData',x,'YData',z);
%%%%%%%%%%%%%%%%%%%%%%%%
section_al = zeros(100,100,100) ; 
for i=1:y_samples
    section=I_use(:,:,i); % get the section data of z_real*x_real at y_real=i
%     figure
%     imshow(x,y,section);
%     axis([0 600 0 700])
    section_act =imresize(section,[100,100]);
 %   imwrite(section_act,['C:\SUHESH\Experiment\PAM\Hair_n_Wire_21st_JUN_13\test2_01\', num2str(i) '.bmp']);
    section_al(i,:,:) = section_act ;
    %imwrite(section_act,['D:\Matlab_files\Alfazal\test2_01_alfazal\', num2str(i) '.bmp']);
    
end
section_am = zeros(50,50,50) ;
for i = 1:2:99
  sect = (section_al(:,:,i)+section_al(:,:,i+1))/2 ;
    sect = imresize(sect,[50,50]) ;
    section_am(:,:,(i+1)/2) = sect ;
end
section_an = zeros(25,25,25) ;
for i = 1:2:49
    sect = (section_am(:,:,i)+section_am(:,:,i+1))/2 ;
    sect = imresize(sect,[25,25]) ;
    section_an(:,:,(i+1)/2) = sect ;
end
section_ap = zeros(12,12,12) ;
for i = 1:2:23
    sect = (section_an(:,:,i)+section_an(:,:,i+1))/2 ;
    sect = imresize(sect,[12,12]) ;
    section_ap(:,:,(i+1)/2) = sect ;
end
    
newsect = section_an ;
for i = 1:25
    for j = 1:25
        for k = 1:25
            newsect(i,j,k) = section_an(j,i,k) ;
            
        end
    end
end




%% Grid Creating
PML_size = 10;  
dx=5*1e-4 ; %[m]
dy=dx ; %[m]
dz=dx ; %[m]
Nx=round(70 *1e-3/dx) ;
Ny=round(70 *1e-3/dy) ;
Nz=round(70 *1e-3/dz) ;
kgrid = makeGrid(Nx, dx, Ny, dy, Nz, dz) ;
medium.sound_speed = 1500; %[m/s]
%Creating time array in k-grid
[kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed);
%% Transducer Details (Inputs needed)
Trans_radius = 36 *1e-3 /dx; %(m)  
Trans_angle = 140 ; %[Degree]
Trans_half_angle = Trans_angle/2 ; %[Degree]
% radius of circle with no sensor on the center of the transducer
radius_cntr_circle = round(5 *1e-3 /dx);
half_angle_missing = (radius_cntr_circle/Trans_radius )*(180/pi) ;
%% Creation of sensor mask (Input needed)
% sphere_offset,x coordinate of where the sensor mask is to be placed 
sphere_offset = 1 ;
Trans_height = round(Trans_radius-(Trans_radius*cosd(70))) ; %[m]
Trans_mask = makeSphericalSection(Trans_radius,Trans_height) ;
sensor.mask = zeros(Nx, Ny, Nz) ;
[a,b,c] = size(Trans_mask) ;

sensor.mask(sphere_offset+(1:a),(Ny-b+1)/2+(1:b),(Nz-c+1)/2+(1:c)) = Trans_mask ;
% x_cntr_msng,y_cntr_msng & z_cntr_msng are the y&z coordinates of the circle where sensor is missing 



   
source.p0 = zeros(Nx,Ny,Nz) ;

%% Defining Source (inputs needed)
source.p0(63:87,60:84,60:84) = newsect ;


max_sr = max(max(max(source.p0)));
for i = 1:Nx
    for j = 1:Ny
        for k=1:Nz
    source.p0(i,j,k) = max_sr/100 ;
        end
    end
end


source.p0 = source.p0/max_sr ;
max_nw = max(max(max(newsect)));
for i = 1:25
    for j = 1:25
        for k = 1:25
            if newsect(i, j, k)< 0.5*max_nw
                newsect(i,j,k) = 0 ;
                
            else
                newsect(i,j,k) = 1 ;
            end
            
            
        end
    end
end


source.p0 = source.p0 * max_sr ;
%%  run the simulation

input_args = {'PMLSize', 10, 'PlotPML', false,'PlotSim', false} ;

sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:}) ;


%%%%%%%%%%%%%  sensor data grouping  %%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                                                                 

%% Finding out the coordinates of the sensor_mask

[x, y, z] = meshgrid( 1:Nx, 1:Ny, 1:Nz ) ;
sensor_mask_cord = [x(:), y(:),z(:), sensor.mask(:)] ;
ind = sensor_mask_cord(:,4)==0 ;
sensor_mask_cord(ind,:) = [] ;
sensor_mask_cord(:,4) = [] ;
sensor_mask_cord(:,[1,2])=sensor_mask_cord(:,[2,1]) ;
sensor_length = length(sensor_mask_cord) ;

%% division of the sensor surface wrt to theta and phi
N_ring = 13 ;
N_snsr = zeros(N_ring,1) ; 
N_snsr(1) = 13 ;
N_snsr(2) = 18 ;
N_snsr(3) = 23 ;
N_snsr(4) = 28 ;
N_snsr(5) = 33 ;
N_snsr(6) = 37 ;
N_snsr(7) = 41 ;
N_snsr(8) = 45 ;
N_snsr(9) = 49 ;
N_snsr(10) = 52 ;
N_snsr(11) = 55 ;
N_snsr(12) = 58 ;
N_snsr(13) = 60 ;



Ntheta_snsr = zeros(N_ring,1) ; 

for i=1:N_ring
    Ntheta_snsr(i) = N_snsr(i) + 1 ;
end

Max_no_theta = max(Ntheta_snsr) ; 
thetas = zeros(1,Max_no_theta) ;
for i = 1:N_ring
    thetas(i,1:Ntheta_snsr(i)) = linspace(0,360,Ntheta_snsr(i)) ;
end
N_snsr_cum = cumsum(N_snsr) ; 
N_theta_cum = cumsum(Ntheta_snsr) ;
thetas_vec = thetas' ;
thetas_vec = thetas_vec(:)' ;
thetas_vec_2 = thetas_vec(thetas_vec~=0) ;

for i = 1:N_ring
    thetas_vec_2 = [thetas_vec_2(1:N_theta_cum(i)-1) 0 thetas_vec_2(N_theta_cum(i):length(thetas_vec_2))] ;
end
thetas_vec_fyn = [0 thetas_vec_2] ;
%thetas_vec_fyn(length(thetas_vec_fyn))  = [] ; 

%%%for i = 1:N_ring
    %%eval(['theta_' num2str(i) '= linspace(0,360,Ntheta_snsr(i))']) ; 
%%end

N_phi = N_ring + 1 ;

% del_theta  = 360/N_theta;
% del_phi  = Trans_half_angle/N_phi;

phi_ele = linspace(half_angle_missing, Trans_half_angle, N_phi) ;
phi_elements = phi_ele.' ;

% total no. of sensor elements
N_angles_total = length(thetas_vec_fyn) ; 





% to specify the sensor element with index no. and the corresponding coordinate in theta and phi
sensor_ele = zeros(N_angles_total, 3);

sensor_ele(:, 1) = 1:N_angles_total ;                                                    % the 1st column is the sensor element index


for i = 1:N_angles_total
    sensor_ele(1:N_angles_total, 2) = thetas_vec_fyn ;       % second column is theta
end


for i = 2:N_ring
    sensor_ele(N_theta_cum(i-1)+1:N_theta_cum(i),3) = phi_elements(i) ;
    sensor_ele(1:N_theta_cum(1),3) = phi_elements(1) ;
end


phi_diff = sensor_ele(16,3)- sensor_ele(2,3) ; 


%% centre of the sphere corresponding to the spherical sensor surface
X0_sens = sphere_offset + Trans_radius + 1 ;
Y0_sens = (Ny-b+1)/2+(b-1)/2+1 ;
Z0_sens = (Nz-c+1)/2+(c-1)/2+1 ;


%% to estimate the polar coordinate of the points on the sensor mask
sensor_mask_cord_R = sqrt((sensor_mask_cord(:, 1) - X0_sens).^2 + (sensor_mask_cord(:, 2) - Y0_sens).^2 + (sensor_mask_cord(:, 3) - Z0_sens).^2 ) ;
sensor_mask_cord_phi = asind(sqrt((sensor_mask_cord(:, 2) - Y0_sens).^2 + (sensor_mask_cord(:, 3) - Z0_sens).^2)./sensor_mask_cord_R) ;  % .*(sensor_mask_cord(:,3)-Z0_sens)./abs(sensor_mask_cord(:,3)-Z0_sens); %for phi negative no need
sensor_mask_cord_theta = atand((( Y0_sens - sensor_mask_cord(:, 2))./((sensor_mask_cord(:, 3) - Z0_sens)))) ;


%% Correction in theta and phi (for 90+ angle in theta and '-'ve angle for phi)


for r = 1:sensor_length;
    if sensor_mask_cord(r,2) <= Y0_sens && sensor_mask_cord(r,3) < Z0_sens
        sensor_mask_cord_theta(r) = sensor_mask_cord_theta(r) + 180 ;
    end
    
    if sensor_mask_cord(r,2) > Y0_sens && sensor_mask_cord(r,3) <= Z0_sens
        sensor_mask_cord_theta(r) = sensor_mask_cord_theta(r) + 180 ;
    end
    
    if sensor_mask_cord(r,2) >= Y0_sens && sensor_mask_cord(r,3) > Z0_sens
        sensor_mask_cord_theta(r) = sensor_mask_cord_theta(r) + 360 ;
    end
    
end

% figure
% plot(sensor_mask_cord_R);
% xlabel('Index No. of Grid Point')
% ylabel('Radius of Grid Points')
% title('Radius of Grid Points')

% figure
% plot(sensor_mask_cord_phi);
% xlabel('Index No. of Grid Point')
% ylabel('Phi Value of Grid Points')
% title('Phi Value of Grid Points')

% figure
% plot(sensor_mask_cord_theta);
% xlabel('Index No. of Grid Point')
% ylabel('Theta Value of Grid Points')
% title('Theta Value of Grid Points')


%% checking of grid points (lying in a particular sensor element) in grouping of grid-points into sensor elements
[m_sensor ,n_sensor] = size(sensor_data) ;
sumll = zeros(N_angles_total, n_sensor+2) ;
count = zeros(N_angles_total, n_sensor+2) ;
tpsensor_data = [sensor_mask_cord_theta sensor_mask_cord_phi sensor_data] ; 

thetas_vec_3 = zeros(N_angles_total,3) ;
thetas_vec_3(1:N_angles_total,:) = sensor_ele ;

%%%%thetas_diff not necessary
thetas_diff = zeros(N_ring,1) ;
for i = 1:N_ring
    %thetas_diff(i) = thetas_vec_3(N_snsr_cum(i),2) -  thetas_vec_3(N_snsr_cum(i)-1,2) ;
end



for j = 1:N_angles_total-1
    tpsensor_data = [sensor_mask_cord_theta sensor_mask_cord_phi sensor_data] ;
    for i = 1:sensor_length
        if (thetas_vec_3 (j, 2) > tpsensor_data(i,1) || tpsensor_data(i,1) >= thetas_vec_3 (j+1, 2))
            tpsensor_data(i,:) = 0 ;
        end
        if (thetas_vec_3(j, 3) > tpsensor_data(i,2)  || tpsensor_data(i,2) >= thetas_vec_3(j, 3)+ phi_diff )
            tpsensor_data(i,:) = 0 ;
        end
    end
    sumll(j,:)= sum(tpsensor_data) ;
    idx = tpsensor_data~=0 ; 
    count(j,:) = sum(idx,1);
end



%% clearing the columns which has the theta and phi values

sumll(:,[1,2]) = [] ;

count(:,[1,2]) = [] ;

%% Average Sensor data

sensr_data_avg = sumll./count ;

%% Making of new sensor mask (Angles) for reconstruction                                                                        

sensor_ele_cntr_angle = zeros(N_angles_total, 3) ;

for j = 1:N_angles_total-1;
    sensor_ele_cntr_angle(j,2) = (sensor_ele(j,2) + sensor_ele(j+1,2))/2 ;
    sensor_ele_cntr_angle(j,3) = sensor_ele(j,3) + (phi_diff)/2 ; 
end

% clearing out the unwanted angles  (rec - rectification)
 sensor_ele_cntr_angle(:, 1) = 1 ;  

for i = 1: N_ring
    sensor_ele_cntr_angle(N_theta_cum(i),:) = 0 ;
end
ind_2 = sensor_ele_cntr_angle(:,1)==0 ;
sensor_ele_cntr_angle(ind_2,:) = [] ;
sensor_ele_cntr_angle(length(sensor_ele_cntr_angle),:) = [] ;
N_snsr_total = length(sensor_ele_cntr_angle) ;

%% Clearing Unwanted Sensor Datas
sensr_data_avg(isnan(sensr_data_avg))=0; 
ind_3 = sensr_data_avg(:,1)==0 ;
sensr_data_avg(ind_3,:) = [] ;


%% Making coordinates of the center of 'sensor elements' from the corresponding angles (Coordinates)
sensor_ele_cntr_cord = zeros(N_snsr_total, 3) ; 
for j = 1:N_snsr_total
    sensor_ele_cntr_cord(j,1) = X0_sens - (Trans_radius * cosd(sensor_ele_cntr_angle(j,3))) ;
    sensor_ele_cntr_cord(j,2) = Y0_sens - (Trans_radius * sind(sensor_ele_cntr_angle(j,3)) * sind(sensor_ele_cntr_angle(j,2))) ;
    sensor_ele_cntr_cord(j,3) = Z0_sens + (Trans_radius * sind(sensor_ele_cntr_angle(j,3)) * cosd(sensor_ele_cntr_angle(j,2))) ;
end
sensor_ele_cntrcord_round = round(sensor_ele_cntr_cord) ;



%% creating binary sensor mask from sensor mask coordinates

sensor.mask = zeros(Nx, Ny, Nz);
for j = 1:N_snsr_total
    sensor.mask(sensor_ele_cntrcord_round(j),sensor_ele_cntrcord_round(j + N_snsr_total),sensor_ele_cntrcord_round(j + 2*N_snsr_total)) = 1 ;
end


%% checking coordinates of sensor_mask

[x, y, z] = meshgrid( 1:Nx, 1:Ny, 1:Nz ) ;
sensor_mask_cord_02 = [x(:), y(:),z(:), sensor.mask(:)] ;
ind = find(sensor_mask_cord_02(:,4)==0) ;
sensor_mask_cord_02(ind,:) = [] ;
sensor_mask_cord_02(:,4) = [] ;
sensor_mask_cord_02(:,[1,2])=sensor_mask_cord_02(:,[2,1]);




%% Arranging sensor data for reconstruction [sensor_data_02] according to the checked coordinates

mask_02_length = length(sensor_mask_cord_02) ;
sensor_data_02 = zeros(mask_02_length ,n_sensor) ;
for i= 1:mask_02_length
    [~,indf]=ismember(sensor_ele_cntrcord_round(i,:),sensor_mask_cord_02,'rows') ;
    sensor_data_02(indf,:) = sensr_data_avg(i,:) ;
end




%% Reconstruction using the new sensor mask and arranged data




kgrid_recon = makeGrid(Nx, dx, Ny, dy, Nz, dz);
kgrid_recon.t_array = kgrid.t_array;
sensor.time_reversal_boundary_data = sensor_data_02;
source.p0 = 0 ;

p0_recon = kspaceFirstOrder3D(kgrid_recon, medium, source, sensor, input_args{:});
%%%optional to load all variables otherwise load ball_1 ball_2






%% to threshold the reconstructed images
p0_recon_mod = p0_recon;
for i = 1:Nx
    for j = 1:Ny
        for k = 1:Nz
            if p0_recon_mod(i, j, k)<0              
                p0_recon_mod(i, j, k) = abs(p0_recon_mod(i, j, k)) ;
            end
        end
    end
end
[m_rec, n_rec, l_rec] = size(p0_recon);
% Fr_no_Sel = 40;
Fr_no_Sel = l_rec/2;
Ln_no_Sel = l_rec/2;

size(p0_recon)
figure,
subplot(2, 2, 1)
imagesc(source.p0(:, :, Fr_no_Sel))
axis('equal')
title('Initial Source')
colorbar 
subplot(2, 2, 3)
imagesc(p0_recon(:, :, Fr_no_Sel))
axis('equal')
title('Reconstructed Image (without threshold)')
colorbar  
subplot(2, 2, 4)
imagesc(p0_recon_mod(:, :, Fr_no_Sel))
title('Reconstructed Image (with threshold)')
axis('equal')
colorbar    


savefig('150source.fig')

for i = 1:Nz 
    imwrite(p0_recon(:, :, i), ['D:\Matlab_files\rslt_25_tra\img_' num2str(i) '.bmp']);
end




for i = 1:l_rec
    imwrite(setup(:, :, i), ['C:\Users\IISERTVM\Documents\MATLAB\New Folder\img_' num2str(i) '.bmp']);
end



%%line ploting & frame by frame thresholding

figure
plot(setup(Ln_no_Sel, :, Fr_no_Sel))
hold on

plot(p0_recon(Ln_no_Sel, :, Fr_no_Sel), '--r')
plot(p0_recon_mod(Ln_no_Sel, :, Fr_no_Sel), ':k')
legend('Original', 'Without Threshold', 'With Threshold')


figure
plot(setup(:, Ln_no_Sel, Fr_no_Sel))
hold on

plot(p0_recon(:, Ln_no_Sel, Fr_no_Sel), '--r')
plot(p0_recon_mod(:, Ln_no_Sel, Fr_no_Sel), ':k')
legend('Original', 'Without Threshold', 'With Threshold')

%To find the maxima on each plane
maxi = zeros(l_rec,1) ;
for i = l_rec
    maxi =  max(max(p0_recon_mod)) ;
end
%Threhold value
threshold_val = 0 ;

thres = maxi * threshold_val ;

p0_thres = p0_recon_mod ;

for i = 1:Nx
    for j = 1:Ny
        for k = 1:Nz
            if p0_thres(i, j, k) < thres(k)
                p0_thres(i, j, k) = 0;
            end
        end
    end
end


