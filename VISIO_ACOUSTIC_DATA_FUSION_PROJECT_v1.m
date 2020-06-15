% ------------------------------------------------------------------------
% VISIO-ACOUSTIC DATA FUSION FOR STRUCTURAL HEALTH MONITORING APPLICATIONS 
% ------------------------------------------------------------------------
%
% ORIGINAL CODE, August 2019
% Original code written by: Chad Samuelson 
%                           Caitrin Duffy-Deno 
%                           Christopher Whitworth
% Reviewed by:  Alessandro Cattaneo 
%               Jeff Tippman 
%               David Mascarenas
%
% INPUT:
%       Hi-FPS video data
%       Hi-Fs microphone data (from 1 microphone)
% OUTPUT:
%       Spectrogram of frequencies in video
%       Location of damaged object in scene emitting hi-frequency sound
%       Modal shape & coordinate plots (if entire damaged object/component is in scene)
%       Image of scene with outlines around each area per frequency
%       Magnified video of object/component of concern
% ------------------------------------------------------------------------
clc
clear all
close all

%% Load data
% Video file
[vid_name,vid_folder]=uigetfile(fullfile(pwd,'folder1','*.mov;*.avi'),...
    'Select a VIDEO file');
if vid_name==0 
    disp('No file selected, exiting...'); 
    return; 
end
v = VideoReader(vid_name);

% Microphone file
[mic_name,mic_folder]=uigetfile(fullfile(pwd,'folder1','*.mat'),...
    'Select a MICROPHONE file');
if mic_name==0 
    disp('No file selected, exiting...'); 
    return; 
end
mic = load(mic_name);
fprintf('Microphone data loaded\n');

framerate=v.FrameRate*v.Duration/.5;
v.CurrentTime=v.Duration/2;
h=waitbar(0,'Loading Video');
nFrames = 1;
while hasFrame(v)
    frame = readFrame(v);
    vid(:,:,:,nFrames) = frame; 
    waitbar(v.CurrentTime/v.Duration,h)
    nFrames = nFrames+1;
    clear frame;
end
nFrames = nFrames-1;
close(h);
h=waitbar(0,'Downsampling Video by 1/5');
nFramesDS=1;
DS_factor = 5;
for ii=1:DS_factor:nFrames
    vidDS(:,:,:,nFramesDS) = vid(:,:,:,ii); 
    waitbar(ii/nFrames,h)
    nFramesDS=nFramesDS+1;
end
nFramesDS=nFramesDS-1;
close(h);
fprintf('Video data loaded and downsampled\n');

%% Constant variable definitions
% Full Frame Sensor PARAMETER
full_fram_diag = 43.27/1000; %m     %Full frame sensor dimensions: 36mm x 24mm 
                                    %therefore diagonal dimension is:
                                    %sqrt(36^2 + 24^2) = 43.27mm

% INPUT
%Camera Specs - Edgertronic SC2 Camera w/50 mm Nikon lens
sensor_w = 17.53/1000; %m
sensor_h = 11.83/1000; %m
pixel_num_w = v.Width; %pixel number
pixel_num_h = v.Height; %pixel number
crop_factor = full_fram_diag/(sqrt(sensor_w^2+sensor_h^2));
focal_length = (50/1000)*crop_factor; %m

%Experimental Set-up info 
% Input metric
%dist_to_object=1.8542; %m (measured 6'1")
% Input ASCI
dist_inches = 93; %inches
dist_m = dist_inches*.0254; %converting (1 in=.0254 m)
dist_to_object = dist_m; %meters

% OUTPUT
%Horizontal Size(top_view)
angular_w = 2*atan(sensor_w/(2*focal_length)); %rad
field_w = 2*dist_to_object*tan(angular_w/2); %m (Lx)
pix_w = field_w/(pixel_num_h*pixel_num_w); %m/pix
%Vertical Size(side_view)
angular_h = 2*atan(sensor_h/(2*focal_length)); %rad
field_h = 2*dist_to_object*tan(angular_h/2); %m (Ly)
pix_h = field_h/(pixel_num_h*pixel_num_w); %m/pix

%% Filter video data, build 3D & 2D Phase and Amplitude Matrices
fprintf('Filtering videos and building Phase & Amplitude Matrices\n');

ht = maxSCFpyrHt(zeros(pixel_num_h,pixel_num_w));
pyrType = 'octave';
switch pyrType
    case 'octave'
        filters = getFilters([pixel_num_h pixel_num_w], 2.^[0:-1:-ht], 4);
        repString = 'octave';
        fprintf('Using octave bandwidth pyramid\n');   
    otherwise 
        error('Invalid Filter Types');
end
[croppedFilters, filtIDX] = getFilterIDX(filters);

buildLevel = @(im_dft, k) ifft2(ifftshift(croppedFilters{k}.* ...
        im_dft(filtIDX{k,1}, filtIDX{k,2})));

fprintf('Moving video to Fourier domain\n');
vidFFT = zeros(pixel_num_h,pixel_num_w,nFrames,'single');
for k = 1:nFrames
    originalFrame = rgb2ntsc(im2single(vid(:,:,:,k)));
    tVid = imresize(originalFrame(:,:,1), [pixel_num_h pixel_num_w]);
    vidFFT(:,:,k) = single(fftshift(fft2(tVid)));
end
vidDSFFT = zeros(pixel_num_h,pixel_num_w,nFramesDS,'single');
for k = 1:nFramesDS
    originalFrame = rgb2ntsc(im2single(vidDS(:,:,:,k)));
    tVid = imresize(originalFrame(:,:,1), [pixel_num_h pixel_num_w]);
    vidDSFFT(:,:,k) = single(fftshift(fft2(tVid)));
end

level = 3; %1=Hi; 2=vert. filt(-); 3=45 cw(\); 4=horiz. filt(l); 5=135 cw(/)...

% Obtain reference phase/frame data
pyrRef = buildLevel(vidFFT(:,:,1), level);
pyrRefDS = buildLevel(vidDSFFT(:,:,1), level);
%pyrRefPhaseOrig = pyrRef./abs(pyrRef);
phaseRef = angle(pyrRef);        
phaseRefDS = angle(pyrRefDS);

fprintf('Building 3D Amplitude & Phase matrices\n');
Phase = zeros(size(pyrRef,1), size(pyrRef,2) ,nFrames,'single');
Amp = zeros(size(pyrRef,1), size(pyrRef,2) ,nFrames,'single');
for frameIDX = 1:nFrames
    filterResponse = buildLevel(vidFFT(:,:,frameIDX), level);
    pyrCurrent = angle(filterResponse);
    Amp(:,:,frameIDX) = abs(filterResponse);
    Phase(:,:,frameIDX) = single(mod(pi+pyrCurrent-phaseRef,2*pi)-pi);                          
end
PhaseDS = zeros(size(pyrRefDS,1), size(pyrRefDS,2) ,nFramesDS,'single');
AmpDS = zeros(size(pyrRefDS,1), size(pyrRefDS,2) ,nFramesDS,'single');
for frameIDX = 1:nFramesDS
    filterResponse = buildLevel(vidDSFFT(:,:,frameIDX), level);
    pyrCurrentDS = angle(filterResponse);
    AmpDS(:,:,frameIDX) = abs(filterResponse);
    PhaseDS(:,:,frameIDX) = single(mod(pi+pyrCurrentDS-phaseRefDS,2*pi)-pi);                          
end 

% fprintf('Temporal Filtering\n');
% temporalFilter = @FIRWindowBP;
% fl = 500; %Hz
% fh = 520; %Hz
% Phase_filt = temporalFilter(Phase, fl/fs,fh/fs);
fprintf('Done\n\n');

%% Anti-Aliasing Idea
fprintf('Perform Anti-Aliasing Algorithm to DownSampled video data\n');

% Transfer DS Video Phase and Mic data into Fourier domain
%Video plots
PhaseDS_FFT=fft(PhaseDS,[],3); %Make 3rd dimension frequency
deltaFFT = zeros(size(PhaseDS_FFT,1)*size(PhaseDS_FFT,2),size(PhaseDS_FFT,3),'single');
Amp_DS_2D = zeros(size(PhaseDS_FFT,1)*size(PhaseDS_FFT,2),size(PhaseDS_FFT,3),'single');
numpix = 1;
h=waitbar(0,'Converting FFT Phase to 2D');
for ii=1:size(PhaseDS_FFT,1)
    for jj=1:size(PhaseDS_FFT,2)
        deltaFFT(numpix,:)=PhaseDS_FFT(ii,jj,:);
        Amp_DS_2D(numpix,:)=AmpDS(ii,jj,:);
        waitbar(numpix/(pixel_num_h*pixel_num_w),h);
        numpix=numpix+1;
    end
end
numpix=numpix-1;
close(h);
sum_pix_signal=zeros(1,size(deltaFFT,2),'single');
h=waitbar(0,'Converting 2D Phase Matrix to 1D Fourier Signal');
for ii=1:numpix
    sum_pix_signal(1,:)=sum_pix_signal+(Amp_DS_2D(ii,:).^2).*deltaFFT(ii,:);
    waitbar(ii/(pixel_num_h*pixel_num_w),h);
end
close(h);
DSframerate=floor(v.Duration*v.FrameRate/(.5*DS_factor));
vidDS_freq=[0:size(PhaseDS_FFT,3)/2-1]*DSframerate/nFramesDS;

% Original Mic. data high-pass filtered @ 120 Hz
fs=mic.Slot4_ai1.Fs; %Hz
filtmic_data = highpass(mic.Slot4_ai1.data,120,fs);
Mic_FFT=fft(filtmic_data,[],1);
mic_freq=[0:(mic.Slot4_ai1.L)/2-1]*mic.Slot4_ai1.Fs/mic.Slot4_ai1.L;

% Subsample mic data
nyq_freq = floor(DSframerate/2);
subnyq_mic = zeros(length(filtmic_data)/nyq_freq);
L_nyqmic = 1; %number of points in subnyquist mic data
for ii=1:mic.Slot4_ai1.L/nyq_freq:mic.Slot4_ai1.L
    subnyq_mic(L_nyqmic) = filtmic_data(ii,1);
    L_nyqmic=L_nyqmic+1;
end
L_nyqmic=L_nyqmic-1;

% Plot FFT's of DS Video & Original Microphone
figure;
subplot(2,1,1); plot(vidDS_freq,abs(sum_pix_signal(1:length(sum_pix_signal)/2))); title(['Fourier Plot of Down Sampled Video @ ',num2str(DSframerate),' frames per second']);
xlabel('Frequency (Hz)');
subplot(2,1,2); plot(mic_freq,abs(Mic_FFT(1:length(Mic_FFT)/2))); title('Fourier Plot of Original Mic Data');
xlabel('Frequency (Hz)');

fprintf('Done\n\n');

%% Produce k-space plots
fprintf('Producing k-space plots\n');
TD = Phase;
TD_DS = PhaseDS;

% Displacement field represented in the fourier domain
FD=fft(TD,[],3); 
FD_pos_freq=FD(:,:,1:size(FD,3)/2);

% Apply 2-FFT to FD to get KD & rearrange K_x & K_y axes
KD=fftshift(fft2(FD_pos_freq)); %Transform from spatial to wave number domain
KD_amp=abs(KD); %Obtain magnitude of KD 

% Plot K-space diagram at certain frequency

% Declare frequency vector
deltaF = fs/nFrames;
freq=[0:(size(FD_pos_freq,3)-1)]*deltaF;

%%%%%%%%%%%%% Automation Way for finding peak frequencies %%%%%%%%%%%%%
% Every pixel will contain the same max frequency...so I chose a random one
%[amp_val freq_ind] = maxk(abs(FD_pos_freq(252,430,2:end)),4,3);
%freq_val = freq_frame*framerate/nFrames;
%freq_val=freq(freq_ind(1,1,2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%% Manual Way for setting peak frequencies %%%%%%%%%%%%%%%%%
freq_val = 510; %Hz
[val idx] = min(abs(freq-freq_val));
freq_ind = idx;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Build Kx axes
DeltaKx=1/field_w;
Kx=[-(pixel_num_w/2):(pixel_num_w/2-1)]*DeltaKx;
%Build Ky axes
DeltaKy=1/field_h;
Ky=[-(pixel_num_h/2):(pixel_num_h/2-1)]*DeltaKy;

%K-space plot
figure;
imagesc(Kx,Ky,KD_amp(:,:,freq_ind)); 
xlabel('K_x','FontSize',24); ylabel('K_y','FontSize',24); title(['K-space Plot @ ',num2str(freq_val),' Hz'],'FontSize',24);
colorbar;
%caxis([.5 2.5]*10^6);
hold on;

%Plot k circle
c=343; %m/s
theta = 0:(pi/36):(2*pi); %rad
k=2*pi*freq_val/c; %rad/m
circ = [k*cos(theta); k*sin(theta)];
plot(circ(1,:),circ(2,:),'r-','LineWidth',2); hold off;






