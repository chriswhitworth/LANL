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

%% Filter video data, build 3D & 2D Phase and Amplitude Matrices
fprintf('Filtering videos and building Phase & Amplitude Matrices\n');
gaborBank = gabor(25,[0 90]); % Combining vertical & horizontal directions

% HI
% Phase=zeros(size(vid,1),size(vid,2),size(vid,4));
% Amp=zeros(size(vid,1),size(vid,2),size(vid,4));
% % Phase_2D=zeros(size(vid,1)*size(vid,2),size(vid,4));
% % Amp_2D=zeros(size(vid,1)*size(vid,2),size(vid,4));
% h=waitbar(0,'Filtering video and building Phase & Amplitude Matrices');
% for ii=1:nFrames
%     vid_ntsc=rgb2ntsc(vid(:,:,:,ii));
%     complex_filt_image = conv2(vid_ntsc(:,:,1),gaborBank(2).SpatialKernel+gaborBank(1).SpatialKernel,'same');
%     Phase(:,:,ii)=angle(complex_filt_image);
%     Amp(:,:,ii)=abs(complex_filt_image);
%     waitbar(ii/nFrames,h);
% end
% close(h);

% h=waitbar(0,'Converting 3D Phase & Amplitude Matrices to 2D');
% numpixels=1;
% for ii=1:size(Phase,1)
%     for jj=1:size(Phase,2)
%         Phase_2D(numpixels,:)=Phase(ii,jj,:);
%         Amp_2D(numpixels,:)=Amp(ii,jj,:);
%         numpixels=numpixels+1;
%         waitbar(numpixels/(size(vid,1)*size(vid,2)),h);
%     end
% end
% numpixels=numpixels-1;
% close(h)

% LO (DOWNSAMPLED)
PhaseDS=zeros(size(vidDS,1),size(vidDS,2),size(vidDS,4));
AmpDS=zeros(size(vidDS,1),size(vidDS,2),size(vidDS,4));
% PhaseDS_2D=zeros(size(vidDS,1)*size(vidDS,2),size(vidDS,4));
% AmpDS_2D=zeros(size(vidDS,1)*size(vidDS,2),size(vidDS,4));
h=waitbar(0,'Filtering downsampled-video and building Phase & Amplitude Matrices');
for ii=1:nFramesDS
    vidDS_ntsc=rgb2ntsc(vidDS(:,:,:,ii));
    complex_filt_image = conv2(vidDS_ntsc(:,:,1),gaborBank(2).SpatialKernel+gaborBank(1).SpatialKernel,'same');
    PhaseDS(:,:,ii)=angle(complex_filt_image);
    AmpDS(:,:,ii)=abs(complex_filt_image);
    waitbar(ii/nFramesDS,h);
end
close(h)

% h=waitbar(0,'Converting 3D Phase & Amplitude Matrices to 2D');
% numpixelsDS=1;
% for ii=1:size(PhaseDS,1)
%     for jj=1:size(PhaseDS,2)
%         PhaseDS_2D(numpixelsDS,:)=PhaseDS(ii,jj,:);
%         AmpDS_2D(numpixelsDS,:)=AmpDS(ii,jj,:);
%         numpixelsDS=numpixelsDS+1;
%         waitbar(numpixelsDS/(size(vidDS,1)*size(vidDS,2)),h);
%     end
% end
% close(h)
% numpixelsDS=numpixelsDS-1;
fprintf('Done\n');

%% Anti-Aliasing Idea
% Transfer Phase and Mic data into Fourier domain
PhaseDS_FFT=fftshift(fft(PhaseDS,[],3),3); %Make 3rd dimension frequency
DSframerate=v.Duration*v.FrameRate*2/DS_factor;
vidDS_freq=[0:size(PhaseDS_FFT,3)-1]*DSframerate/nFramesDS;

Mic_FFT=fftshift(fft(mic.Slot5_ai2.data,[],1));
mic_freq=[0:mic.Slot5_ai2.L-1]*mic.Slot5_ai2.Fs/mic.Slot5_ai2.L;

% Plot FFT's of DS Video & Microphone





