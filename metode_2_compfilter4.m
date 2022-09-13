%% Metode: Complementary Filter (CF)

%% Komputasi Respiration Rate Dengan CF adalah metode
%% fusi antara Euler angle yang dihasilkan accelerometer
%% dan gyroscope angle. 
%% Data mentah Acc dan Gyr difilter terlebih dahulu dengan
%% BPF 0.1 sd 0.8 Hz untuk mendapatkan frequensi respirasi.

%% Respiration rate dideteksi berdasarkan angular displacement dari gerakan
%% pitch chest wall surface.
%% Sensor IMU diletakan di sekitar upper abdomen - sternum
%% sensor pada upper abdomen menghasilkan sinyal yang
%% bagus dibandingkan dengan peletakan di sternum.
%% Berikut algoritmanya:

%%  1a.Acc --> 2.BPS --> 3.Euler --> 4.Complementary --> 5.MVA --> 6.Peak Detection
%%                                            ^
%%                                            |  
%%  1b.Gyr -----------------------------------

clc;
clear all;
close all;

%% Read Data from excel

%% Dataset 1: Abdomen Supine Position 
%% Utk data abdomen, freq dominan 0.2< f< 0.8 Hz pada FFT sumbu X
%% Utk abdomen, freq dominan <1 Hz pada FFT

dataset='01_abdomen_13bpm.xls';% bpm = 12, brv= 4.77;
%dataset='02_abdomen_13bpm.xls';% bpm = 13, brv= 4.55;
%dataset='03_abdomen_11bpm.xls';% bpm = 11 brv= 5.3;
%dataset='04_abdomen_11bpm.xls';% bpm = 12 brv= 5.07
%dataset='05_abdomen_14bpm.xls';%bpm = 14 brv= 4.3
%dataset='09-abdomen-13rpm.xls'; %rpm = 12, brv= 4.93;
%dataset='08-abdomen-12bpm.xls'; %rpm = 11, brv= 5.16;

%dataset='13_abdomen_12rpm.xls' %rpm = 12 rrv=4.88%
%dataset='14_abdomen_12rpm.xls' %rpm = 11 rrv=5.37
%dataset='16_abdomen_14rpm.xls'  %rpm = 14, rrv=4.33

%dataset='01-RR-13rpm_abdomen.xls'%rpm = 12, rrv=5.01
%dataset='02-RR-12rpm-abdomen.xls'%rpm = 11, rrv=5.37

%dataset='03-RR-13rpm-abdomen.xls'%rpm = 13, rrv=5.37

%% Dataset 1: Sternum Supine Position ------------------------------------
%% Utk Sternum, freq dominan 3< f<15 Hz pada FFT

%dataset='01_sternum_12bpm.xls';%rpm=12 rrv= 4,95
%dataset='04_sternum13bpm.xls'; %rpm=13 rrv= 4,48
%dataset='05_sternum13bpm.xls'; %rpm=11 rrv= 5,18

%dataset='06_sternum_13rpm.xls' %rpm=11 rrv= 5,18 ; roll
%dataset='07_sternum_12rpm.xls'%rpm=12 rrv= 4,72
%dataset='08_sternum_12rpm.xls' %%rpm=11 rrv= 5,016 ; roll


%dataset='17_sternum_13rpm.xls'% rpm=12 rrv= 4,84 ; roll
%dataset='11_sternum_11bpm.xls' %% rpm=9 rrv= 6,53 roll
%dataset='12_sternum_11bpm.xls' % rpm = 10, brv= 4.13

%dataset='01-RR1-13rpm_sternum.xls'%rpm = 11, brv= 5.07 roll
%dataset='02-RR-12rpm-sternum.xls'% rpm = 11, brv= 4.13 roll

%dataset='03-RR-13rpm-sup.xls'% rpm = 11, brv= 4.13 roll

%----------- Pembacaan Data Dari Sensor MBient IMU -----------%
data = xlsread(dataset);
disp('File Name:'); 
disp(dataset);

crop=1:6000; % crop signal to remove unintended initial noise
time=data(crop,1);
%t = time;

%% read accelerometer data
accx=data(crop,6);%*9.81; %convert g to /ms2, if data in g
ax = accx;

accy=data(crop,7);%*9.81;% acc y
ay = accy;

% khusus acc z
accz=data(crop,8);%*9.81;%convert g to /ms2, if data in g
az = accz;
az = az-9.81 ; % Posisi utk telentang, az hrs dikurangi gravitasi

% baseline removal (DC Component removal)
ax=ax- (sum(ax)./size(ax,1));
ay=ay-(sum(ay)./size(ay,1));
az=az- (sum(az)./size(az,1));

%% read Gyroscope data
gx=data(crop,9);% gyr x
gy=data(crop,10);% gyr y
gz=data(crop,11);% gyr z

fs =1/(data(5,1)-data(4,1));
t= (0:length(time)-1)/fs;
dt=data(5,1)-data(4,1);
batasx=[0 60];
%% BPF
%% BPF on Accelerometer Az axis
  [b, a] = butter(1, [0.1/fs 0.8/fs], 'bandpass');
  axBPF = filtfilt(b, a, ax );

%% BPF on Accelerometer Az axis
  [b, a] = butter(1, [0.1/fs 0.8/fs], 'bandpass');
  ayBPF = filtfilt(b, a, ay );
  
%% BPF on Accelerometer Az axis
  [b, a] = butter(1, [0.1/fs 0.8/fs], 'bandpass');
  azBPF = filtfilt(b, a, az );
  
%% Convert gyroscope measurements to radians

Gx_rad = gx * pi / 180.0;
Gy_rad = gy * pi / 180.0;
Gz_rad = gz * pi / 180.0;
   

%% %% BPF Based Complementary
%% 1) BPF-Based Complementary Filter: EULER Angle formula
%% a) Obtain Euler angle of pitch and roll from Accelerometer
tic; % start recording time
phi_hat_acc   = atan2(ayBPF, sqrt(axBPF .^ 2 + azBPF .^ 2)); % if roll is used
theta_hat_acc = atan2(-axBPF, sqrt(ayBPF .^ 2 + azBPF .^ 2));% if pitch is used

%% b) integral of gyroscope (degree/s) to obtain gyro angle (degree)
phi_hat_gyr   = zeros(1, length(t)); % roll
theta_hat_gyr = zeros(1, length(t)); % pitch

for i = 2:length(t)
   p = Gx_rad(i);
   q = Gy_rad(i);
   r = Gz_rad(i);
   
   phi_hat   = phi_hat_gyr(i - 1);
   theta_hat = theta_hat_gyr(i - 1);
    
   phi_hat_gyr(i)   = phi_hat   + dt * (p + sin(phi_hat) * tan(theta_hat) * q + cos(phi_hat) * tan(theta_hat) * r);
   theta_hat_gyr(i) = theta_hat + dt * (cos(phi_hat) * q - sin(phi_hat) * r);
   
end

%% 2) Complimentary Filter from Euler angle

alpha = 0.2; % constant of alpha

phi_hat_complimentary   = zeros(1, length(t));
theta_hat_complimentary = zeros(1, length(t));

for i=2:length(t)
    p = Gx_rad(i);
    q = Gy_rad(i);
    r = Gz_rad(i);
   
    phi_hat   = phi_hat_complimentary(i - 1);
    theta_hat = theta_hat_complimentary(i - 1);
    
    phi_hat_gyr_comp   = phi_hat   + dt * (p + sin(phi_hat) * tan(theta_hat) * q + cos(phi_hat) * tan(theta_hat) * r);
    theta_hat_gyr_comp = theta_hat + dt * (cos(phi_hat) * q - sin(phi_hat) * r);
       
    phi_hat_complimentary(i)   = (1 - alpha) * phi_hat_gyr_comp   + alpha * phi_hat_acc(i);
    theta_hat_complimentary(i) = (1 - alpha) * theta_hat_gyr_comp + alpha * theta_hat_acc(i);    
end

elapsed2=toc;                              % stop recording time
disp('Complementary Filter Computation:'); % to display in command line
disp(elapsed2);                            % to display in command line
%% ---------------------------------------

%% 4. Remove signal drift of angular value
filtCutOff = 0.5;
[b, a] = butter(1, (filtCutOff)/(fs), 'high');
theta_hat_complimentary = filtfilt(b, a, theta_hat_complimentary);
phi_hat_complimentary = filtfilt(b, a, phi_hat_complimentary);

%% END OF COMPUTATION

 
snr_phi=round(snr(phi_hat_complimentary),2) ;
snr_theta=round(snr(theta_hat_complimentary),2);

if snr_theta >= snr_phi
        fusion= theta_hat_complimentary;
        sudut = '  Pitch';
    else
        fusion= phi_hat_complimentary;
        sudut = '  Roll';
end


%% Figure 5: Results BPF based Complementary Filter
figure('Position', [800 500 800 600], 'NumberTitle', 'off', 'Name', 'BPF-Metode Complementary Filter');    
ax(1) = subplot(2,1,1);
    hold on;
    
    
    tic
   
    
    %% BPF on Accelerometer Az axis
    [b, a] = butter(1, [0.2/fs 0.8/fs], 'bandpass');
    out = filtfilt(b, a, fusion );
    
    elapsed3=toc;
    disp('Moving average Computation:'); 
    disp(elapsed3);
    
    
    %plot(time, pitchAngle, 'r');hold on;
    %plot(time, out, 'g','Linewidth',4);hold on;
    plot(time, out, 'b','Linewidth',2);hold on;
    
    %title('Respiration Signal derived from ', sudut);
    
    title({['Respiration Signal derived from',sudut, ' Angle']
           [' Dataset : ',dataset]  
            },'Interpreter','none');
    
    xlabel('Time (s)');
    ylabel('Degree (^{o})');
    xlim(batasx);
    %ylim([-0.2 0.2]);
    hold off;
    set(gcf,'color','w');
    set(gca,'FontSize',16);grid on;
    ax = gca; 
    ax.GridLineStyle = '--'; 
    ax.GridColor = 'k'; 
    ax.GridAlpha = 0.25;
    ax.LineWidth=1;
    ax.Color=[0.90 0.90 0.90];
    set(gcf,'color','w');
    grid minor;
    grid on;
    
    
ax(2) = subplot(2,1,2);
   
    %% 4. Baseline removal of complementary filter
    %%    output
    tic
    
    elapsed4=toc;
    disp('Baseline Removal Computation:'); 
    disp(elapsed4);
    
    windowWidth = 75; % jumlah data dalam window.
    kernel = ones(windowWidth,1) / windowWidth;
    rrdata = out; % tkae the moving average filter data of comp filter
    thres = 0.3*max(rrdata); % 30% of max thres
    peakdist=2.75; %antara 3 sd 4 second
    minthres = 0.3*min(rrdata);
    %rrdata(rrdata<0) = 0;
    kernel = ones(windowWidth,1) / windowWidth;
    rrdata = filter(kernel, 1, rrdata);
    
    %% 5. Peak detection
    tic
    [rrpeaks,locs] = findpeaks(rrdata,t,'MinPeakHeight',thres,'MinPeakDistance',peakdist);
    [minpeaks,minlocs] = findpeaks(-rrdata,t,'MinPeakHeight',minthres,'MinPeakDistance',3.5);
    pks = numel(rrpeaks);
    elapsed5=toc;
    disp('Peak Detection Computation:'); 
    disp(elapsed5);
    
    
    hold on;
    
    plot(locs,rrpeaks,'kO','Linewidth',4);hold on
    text(locs+.02,rrpeaks,num2str((1:numel(rrpeaks))'),'FontSize',18, 'FontWeight','bold', 'HorizontalAlignment','left', 'VerticalAlignment','middle')
    plot(t, rrdata, 'b','Linewidth',2);hold on
    %H=area(t,rrdata);
    %set(H,'FaceColor',[1 0.5 0]);hold on;
    
    
    %plot(minlocs,-minpeaks,'k*','Linewidth',3);hold on
    xlabel('Time (s)');
    xlim(batasx);
    %legend('Fiducial Points', 'Filtered Signal','Min point');
    ylabel('Degree (^{o})');
    y1=get(gca,'ylim');
    rrv=round(mean(diff(locs)),2);
    
    title({
           [' Estimated Respiration Rate per Minute : ',num2str(pks),' rpm']
           [' Respiration Rate Variation : ',num2str(rrv),' sec']  
           [' Dataset : ',dataset]
            },'Interpreter','none');
    
    hold off;
    set(gca,'FontSize',16);grid on;
    ax = gca; 
    ax.GridLineStyle = '--'; 
    ax.GridColor = 'k'; 
    ax.GridAlpha = 0.25;
    ax.LineWidth=1;
    ax.Color=[0.90 0.90 0.90];
    set(gcf,'color','w');
    grid minor;
    grid on;  
    