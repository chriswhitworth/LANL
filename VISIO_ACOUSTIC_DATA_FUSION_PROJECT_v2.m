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
% Microphone file
[mic_name,mic_folder]=uigetfile(fullfile(pwd,'folder1','*.mat'),...
    'Select a MICROPHONE file');
if mic_name==0 
    disp('No file selected, exiting...'); 
    return; 
end
mic = load(mic_name);
fprintf('Microphone data loaded\n');

% Video file
[vid_name,vid_folder]=uigetfile(fullfile(pwd,'folder1','*.mov;*.avi'),...
    'Select a VIDEO file');
if vid_name==0 
    disp('No file selected, exiting...'); 
    return; 
end
v = VideoReader(vid_name);

v.CurrentTime=v.Duration/2; %start half way through video
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
fprintf('Video Data Loaded\n');

%% Filter video data, build 3D & 2D Phase and Amplitude Matrices
fprintf('Filtering video and building Phase & Amplitude Matrices\n');

ht = maxSCFpyrHt(zeros(v.Height,v.Width));
pyrType = 'octave';
switch pyrType
    case 'octave'
        filters = getFilters([v.Height v.Width], 2.^[0:-1:-ht], 4);
        repString = 'octave';
        fprintf('Using octave bandwidth pyramid\n');   
    otherwise 
        error('Invalid Filter Types');
end
[croppedFilters, filtIDX] = getFilterIDX(filters);

buildLevel = @(im_dft, k) ifft2(ifftshift(croppedFilters{k}.* ...
        im_dft(filtIDX{k,1}, filtIDX{k,2})));

fprintf('Moving video to Fourier domain\n');
vidFFT = zeros(v.Height,v.Width,nFrames,'single');
for k = 1:nFrames
    originalFrame = rgb2ntsc(im2single(vid(:,:,:,k)));
    tVid = imresize(originalFrame(:,:,1), [v.Height v.Width]);
    vidFFT(:,:,k) = single(fftshift(fft2(tVid)));
end

level = 3; %1=Hi; 2=vert. filt(-); 3=45 cw(\); 4=horiz. filt(l); 5=135 cw(/)...

% Obtain reference phase/frame data
pyrRef = buildLevel(vidFFT(:,:,1), level);
%pyrRefPhaseOrig = pyrRef./abs(pyrRef);
phaseRef = angle(pyrRef);

fprintf('Building 3D Amplitude & Phase matrices\n');
Phase = zeros(size(pyrRef,1), size(pyrRef,2) ,nFrames,'single');
Amp = zeros(size(pyrRef,1), size(pyrRef,2) ,nFrames,'single');
for frameIDX = 1:nFrames
    filterResponse = buildLevel(vidFFT(:,:,frameIDX), level);
    pyrCurrent = angle(filterResponse);
    Amp(:,:,frameIDX) = abs(filterResponse);
    Phase(:,:,frameIDX) = single(mod(pi+pyrCurrent-phaseRef,2*pi)-pi);                          
end

% fprintf('Temporal Filtering\n');
% temporalFilter = @FIRWindowBP;
% fl = 500; %Hz
% fh = 520; %Hz
% Phase_filt = temporalFilter(Phase, fl/fs,fh/fs);
fprintf('Done\n\n');

%% Anti-Aliasing Idea
fprintf('Perform Anti-Aliasing Algorithm to DownSampled video data\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% ANTI-ALIASING STEPS %%%%%%%%%%%%%%%%%%%%%%%%%
% 1. DOWN SAMPLE VIDEO DATA
%       1.a)    Down Sample video data
%       1.b)    Plot Hi-FPS vs Lo-FPS video data fft plots
% 2. DOWN SAMPLE MICROPHONE DATA
%       2.a)    Filter out frequencies below 100 Hz on orig mic data
%       2.b)    Down Sample Mic data
%       2.c)    Plot Hi-Fs vs Lo-Fs mic data fft plots
% 3. OVERLAY AND MATCH FFT PLOTS
%       3.a)    Use findpeaks function to find peaks in mic data
% 4. Make DS Alias calculation & Filter DS video accordingly
%       4.a)    Create alias finding calculation on found peak frequencies
%       4.b)    Filter DS video data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 1.a: DOWN SAMPLE VIDEO DATA
DS_factor = 5; %downsampling video by 1/5
fprintf(['Down Sampling Video Data by 1/',num2str(DS_factor),'\n']);
PhaseDS = zeros(size(Phase,1),size(Phase,2),ceil(size(Phase,3)/DS_factor));
AmpDS = zeros(size(Amp,1),size(Amp,2),ceil(size(Amp,3)/DS_factor));
nFramesDS = 1;
for ii=1:DS_factor:nFrames
    AmpDS(:,:,nFramesDS)=Amp(:,:,ii);
    PhaseDS(:,:,nFramesDS)=Phase(:,:,ii);
    nFramesDS=nFramesDS+1;
end
nFramesDS=nFramesDS-1;
fprintf('Done\n');

%% STEP 1.b: PLOT HIGH-FPS VIDEO FOURIER SPECTRUM VS. DOWN SAMPLED VIDEO DATA
fprintf('Plotting Hi-FPS Video vs Lo-FPS Video spectrums\n');
% HI-FPS VIDEO
Phase_FFT = fft(Phase,[],3); %Convert 3rd dimension of phase to be frequency
fourier_vid = zeros(size(Phase,3),1,'single');
h=waitbar(0,'Converting 3D Phase Matrix to 1D Fourier Signal');
for ii=1:nFrames
    weighted_phase=sum(sum(abs((Amp(:,:,ii)).^2).*Phase_FFT(:,:,ii)));
    sumamp=sum(sum(abs(Amp(:,:,ii))));
    fourier_vid(ii)=mean(weighted_phase)/sumamp;
    waitbar(ii/nFrames,h);
end
close(h);

%LO-FPS VIDEO
PhaseDS_FFT = fft(PhaseDS,[],3); %Convert 3rd dimension of phase to be frequency
fourier_vidDS = zeros(size(PhaseDS,3),1,'single');
h=waitbar(0,'Converting 3D Down Sampled Phase Matrix to 1D Fourier Signal');
for ii=1:nFramesDS
    weighted_phase=sum(sum(abs((AmpDS(:,:,ii)).^2).*PhaseDS_FFT(:,:,ii)));
    sumamp=sum(sum(abs(AmpDS(:,:,ii))));
    fourier_vidDS(ii)=mean(weighted_phase)/sumamp;
    waitbar(ii/nFramesDS,h);
end
close(h);

% PLOT to compare hi-fps to lo-fps video data
% Set-up Frequency Vectors
framerate=v.FrameRate*v.Duration/.5; %.5s=how long video was actually recorded
vid_freq=[0:size(Phase_FFT,3)/2-1]*framerate/nFrames;
vid_pos=abs(fourier_vid(1:length(fourier_vid)/2));
DSframerate=framerate/DS_factor;
vidDS_freq=[0:size(PhaseDS_FFT,3)/2-1]*DSframerate/nFramesDS;
vidDS_pos=abs(fourier_vidDS(1:length(fourier_vidDS)/2));

% Plot
figure;
subplot(2,1,1); 
plot(vid_freq,vid_pos); 
title(['Hi-Frame Per Second Video Data Spectrum @ ',num2str(framerate),' frames per second'],'FontSize',20);
xlabel('Frequency (Hz)','FontSize',20);
subplot(2,1,2);
plot(vidDS_freq,vidDS_pos); 
title(['Lo-Frame Per Second Video Data Spectrum @ ',num2str(DSframerate),' frames per second'],'FontSize',20);
xlabel('Frequency (Hz)','FontSize',20);
fprintf('Done\n\n');

%% STEP 2.a: Filter out frequencies below 100 Hz
micfilt=highpass(mic.Slot4_ai1.data,100,mic.Slot4_ai1.Fs); 

%% STEP 2.b: DOWN SAMPLE MICROPHONE DATA
fprintf('Down Sampling Microphone Data\n');
% Subsample mic data
micDS = zeros(mic.Slot4_ai1.L/floor(mic.Slot4_ai1.Fs/DSframerate),1);
Length_nyqmic = 1; %number of points in subnyquist mic data
for ii=1:floor(mic.Slot4_ai1.Fs/DSframerate):mic.Slot4_ai1.L
    micDS(Length_nyqmic) = micfilt(ii,1);
    Length_nyqmic=Length_nyqmic+1;
end
Length_nyqmic=Length_nyqmic-1;
fprintf('Done\n');

%% STEP 2.c: Plot Hi-Fs vs Lo-Fs mic data fft plots
fprintf('Plotting Hi-FS mic vs Lo-FS mic data\n');
%Transform mic to fourier domain
micfreq = mic.Slot4_ai1.Fs*[0:(mic.Slot4_ai1.L-1)]/(mic.Slot4_ai1.L);
micFFT = fft(micfilt);
micfreq_pos = micfreq(1:length(micfreq)/2);
micFFT_pos = abs(micFFT(1:length(micFFT)/2));
micDSfreq = DSframerate*[0:(Length_nyqmic-1)]/Length_nyqmic;
micDSFFT = fft(micDS);
micDSfreq_pos = micDSfreq(1:length(micDSfreq)/2);
micDSFFT_pos = abs(micDSFFT(1:length(micDSFFT)/2));
% PLOT to compare hi-fs to lo-fs mic data
figure;
subplot(2,1,1); 
plot(micfreq_pos,micFFT_pos); 
title(['Hi-Sampling Rate Mic Data Spectrum @ ',num2str(5000),' Sampling Rate'],'FontSize',20);
xlabel('Frequency (Hz)','FontSize',20);
subplot(2,1,2);
plot(micDSfreq_pos,micDSFFT_pos); 
title(['Lo-Sampling Rate Mic Data Spectrum @ ',num2str(DSframerate),' Sampling Rate'],'FontSize',20);
xlabel('Frequency (Hz)','FontSize',20);
fprintf('Done\n\n');

%% STEP 3: OVERLAY AND MATCH FFT PLOTS
%% STEP 3.a: Use findpeaks function to find peaks in mic data
%% Mic data
% Original Mic data
[PKSmic,location_mic_Hz]=findpeaks(micFFT_pos,micfreq_pos,'MinPeakHeight',.4,...
    'MinPeakDistance',20,'SortStr','ascend','NPeaks',10);
PKSmic;
location_mic_Hz
figure;
subplot(2,1,1);
findpeaks(micFFT_pos,micfreq_pos,'MinPeakHeight',.4,...
    'MinPeakDistance',20,'SortStr','none','NPeaks',10,'Annotate','peaks');
title('Peaks found in original mic data (at least 200 Hz apart)','FontSize',20);
xlabel('Frequency (Hz)','FontSize',20);

% Down Sampled Mic data
% [PKSDSmic,indices_DSmic]=findpeaks(micDSFFT_pos,'MinPeakHeight',.09,...
%     'MinPeakDistance',20,'SortStr','ascend','NPeaks',10);
% indices_DSmic
% [PKSDSmic,location_DSmic_Hz]=findpeaks(micDSFFT_pos,micDSfreq_pos,'MinPeakHeight',.09,...
%     'MinPeakDistance',20,'SortStr','ascend','NPeaks',10);
% PKSDSmic;
% location_DSmic_Hz
% subplot(2,1,2);
% findpeaks(micDSFFT_pos,micDSfreq_pos,'MinPeakHeight',.09,...
%     'MinPeakDistance',20,'SortStr','none','NPeaks',10);
% title('Peaks found in Down Sampled mic data (at least 20 Hz apart)','FontSize',20);
% xlabel('Frequency (Hz)','FontSize',20);

%% STEP 4: Make DS Alias calculation & Filter DS video accordingly
%% STEP 4.a:    Create alias finding calculation on found peak frequencies
alias_freq=zeros(length(location_mic_Hz));
for ii=1:length(location_mic_Hz)
    currFreq=location_mic_Hz(ii);
    currInt=floor(currFreq/(DSframerate/2));
    if mod(currInt,2)==0 %Even
        
    else %Odd
        
    end
end


%% STEP 4.b:    Bandpass-Pass Filter video data based on mic data peaks






