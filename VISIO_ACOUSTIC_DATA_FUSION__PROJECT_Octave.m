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
[vid_name,vid_folder]=uigetfile(fullfile(pwd,'folder1','*.mov'),...
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

v.CurrentTime=0;
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
nFrameDS=1;
DS_factor = 5;
for ii=1:DS_factor:nFrames
    vid_downsample(:,:,:,nFrameDS) = vid(:,:,:,ii); 
    waitbar(ii/nFrames,h)
    nFrameDS=nFrameDS+1;
end
nFrameDS=nFrameDS-1;
close(h);
fprintf('Video data loaded and downsampled\n');

%% Get user input
% Create user input dialogue box
prompt = {'Analysis Direction (x, 45ccw, y, or 135ccw)','Number of Modes to Identify','Reference Frame'};
dlg_title = 'Processing Parameters';
default_params ={'x','3','1'};
params = inputdlg(prompt,dlg_title,1,default_params);

direction = params{1}; %Filter direction
nmode = str2double(params{2}); %Number of modes
refFrameIDX = str2double(params{3}); %Reference frame

%% Build pyramid from video data - MAKE EDITS TO US ENTIRE PYRAMID -

ht = maxSCFpyrHt(zeros(v.Height,v.Width));
pyrType = 'octave'; % Can change this variable
switch pyrType
    case 'octave'
        filters = getFilters([v.Height v.Width], 2.^[0:-1:-ht], 4);
        repString = 'octave';
        fprintf('Using octave bandwidth pyramid\n');        
    case 'halfOctave'            
        filters = getFilters([v.Height v.Width], 2.^[0:-0.5:-ht], 8,'twidth', 0.75);
        repString = 'halfOctave';
        fprintf('Using half octave bandwidth pyramid\n'); 
    case 'smoothHalfOctave'
        filters = getFiltersSmoothWindow([v.Height v.Width], 8, 'filtersPerOctave', 2);           
        repString = 'smoothHalfOctave';
        fprintf('Using half octave pyramid with smooth window.\n');
    case 'quarterOctave'
        filters = getFiltersSmoothWindow([v.Height v.Width], 8, 'filtersPerOctave', 4);
        repString = 'quarterOctave';
        fprintf('Using quarter octave pyramid.\n');
    otherwise 
        error('Invalid Filter Types');
end
[croppedFilters, filtIDX] = getFilterIDX(filters);
numLevels = numel(filters);

% Function converts frame back to physical space from Fourier
%   after applying the spatial filters
buildLevel = @(im_dft,lvl) ifft2(ifftshift(croppedFilters{lvl}.* ...
    im_dft(filtIDX{lvl,1}, filtIDX{lvl,2})));

% Moving to fourier domain and then applying pyramid to each frame
h = waitbar(0,'Moving hi-FPS video into Fourier domain...');
vidFFT = zeros(v.Height,v.Width,nFrames,'single');
for k = 1:nFrames
    ntsc_frame = rgb2ntsc(im2single(vid(:,:,:,k)));
    lim_ntsc_frame = imresize(ntsc_frame(:,:,1),[v.Height v.Width]);
    vidFFT(:,:,k) = single(fftshift(fft2(lim_ntsc_frame)));
    waitbar(k/nFrames,h);   
end
close(h)
h = waitbar(0,'Moving low-FPS video into Fourier domain...');
vidFFT_DS = zeros(v.Height,v.Width,nFrameDS,'single');
for k = 1:nFrameDS
    ntsc_frame = rgb2ntsc(im2single(vid_downsample(:,:,:,k)));
    lim_ntsc_frame = imresize(ntsc_frame(:,:,1),[v.Height v.Width]);
    vidFFTDS(:,:,k) = single(fftshift(fft2(lim_ntsc_frame)));
    waitbar(k/nFrameDS,h);   
end
close(h)

switch direction
    case 'x',       level=2; % X-direction filters
    case '45ccw',   level=3; % 45 degrees counter-clockwise from x-direction
    case 'y',       level=4; % Y-direction filters (90 deg ccw)
    case '135ccw',  level=5; % 135 degrees counter-clockwise from x-direction
    otherwise,  error('Only Enter x, 45ccw, y, or 135ccw for Analysis Direction')
end

% ------------------------------------------------------------------------
%HI-FPS

pyrRef = buildLevel(vidFFT(:,:,refFrameIDX),level);
pyrRef = angle(pyrRef);


% Initialize matrices
Phase = zeros(size(pyrRef,1), size(pyrRef,2) ,nFrames,'single');
Amp = zeros(size(pyrRef,1), size(pyrRef,2) ,nFrames,'single');

% Compute Phase and amplitude for each frame, make phases
%   unwrapped and relative to reference frame 
ha = waitbar(0,'Computing Filtered Phases HI-FPS video...');
for frameIDX = 1:nFrames
    filterResponse = buildLevel(vidFFT(:,:,frameIDX),level);
    pyrCurrent = angle(filterResponse);
    % Local amplitude matrix
    Amp(:,:,frameIDX)= abs(filterResponse); 
    % Local phase matrix (shifted to within [-pi,pi])
    Phase(:,:,frameIDX) = single(mod(pi+pyrCurrent-pyrRef,2*pi)-pi); 
    waitbar(frameIDX/nFrames,ha)
end
close(ha)

% ------------------------------------------------------------------------
%LOW-FPS - DOWNSAMPLED(DS)

pyrRefDS = buildLevel(vidFFTDS(:,:,refFrameIDX),level);
pyrRefDS = angle(pyrRefDS);

% Initialize matrices
PhaseDS = zeros(size(pyrRefDS,1), size(pyrRefDS,2) ,nFrameDS,'single');
AmpDS = zeros(size(pyrRefDS,1), size(pyrRefDS,2) ,nFrameDS,'single');

% Compute Phase and amplitude for each frame, make phases
%   unwrapped and relative to reference frame 
ha = waitbar(0,'Computing Filtered Phases LOW-FPS video...');
for frameIDX = 1:nFrameDS
    filterResponseDS = buildLevel(vidFFTDS(:,:,frameIDX),level);
    pyrCurrentDS = angle(filterResponseDS);
    % Local amplitude matrix
    AmpDS(:,:,frameIDX)= abs(filterResponseDS); 
    % Local phase matrix (shifted to within [-pi,pi])
    PhaseDS(:,:,frameIDX) = single(mod(pi+pyrCurrentDS-pyrRefDS,2*pi)-pi); 
    waitbar(frameIDX/nFrameDS,ha)
end
close(ha)

%% Obtain 2D phase & amplitude whose rows=frame_# and cols=pixel_#(going across then down)

fprintf('Building 2-D matrices from 3-D matrices...\n');

% ------------------------------------------------------------------------
%HI-FPS
Phase_2D = zeros(size(Phase,3),size(Phase,1)*size(Phase,2));
thresh = multithresh(mean(Amp(:,:,1:end),3));

pixel_num=1;
ll=1;
for ii = 1:size(Phase,1)
    for jj = 1:size(Phase,2)
        Phase_iijj=Phase(ii,jj,:);
        
        % Reshape pixels at high intensity to a time series
        %   along the column space of Ph_2D and save coordinate
        %   locations in the frame to coord_list for future
        %   reconstruction.
        if mean(Amp(ii,jj,1:end))>=thresh
            Phase_2D(:,pixel_num)=Phase_iijj;
            coord_list(pixel_num,:)=int32([ii,jj]);
            pixel_num = pixel_num + 1;
        else
            % Put noise into the recon matrix with spaces held 
            %   for signal by columns of ones
            Phase_2D_recon(:,ll)=Phase_iijj;
        end
        ll = ll+1;
    end
end
pixel_num = pixel_num - 1;

% ------------------------------------------------------------------------
%LOW-FPS - DOWNSAMPLED(DS)
Phase_2D_DS = zeros(size(PhaseDS,3),size(PhaseDS,1)*size(PhaseDS,2));
threshDS = multithresh(mean(AmpDS(:,:,1:end),3));

pixel_numDS=1;
ll=1;
for ii = 1:size(PhaseDS,1)
    for jj = 1:size(PhaseDS,2)
        PhaseDS_iijj=PhaseDS(ii,jj,:);
        
        % Reshape pixels at high intensity to a time series
        %   along the column space of Ph_2D and save coordinate
        %   locations in the frame to coord_list for future
        %   reconstruction.
        if mean(Amp(ii,jj,1:end))>=thresh
            Phase_2D_DS(:,pixel_numDS)=PhaseDS_iijj;
            coord_listDS(pixel_numDS,:)=int32([ii,jj]);
            pixel_numDS = pixel_numDS + 1;
        else
            % Put noise into the recon matrix with spaces held 
            %   for signal by columns of ones
            Phase_2D_reconDS(:,ll)=PhaseDS_iijj;
        end
        ll = ll+1;
    end
end
pixel_numDS = pixel_numDS - 1;
fprintf('Done\n');

%% Understanding Non-negative Matrix Factorization

[W,H] = nnmf(Phase_2D,2);
%% Apply PCA/CP to 2D Phase Matrices, then Build Mode Shape and Modal Coordinate matrices

fprintf('Applying PCA to 2D phase matrices...\n');
% ------------------------------------------------------------------------
%HI-FPS
[COEFF, SCORE, ~, ~, explained] = pca(Phase_2D,'Algorithm','svd');
%nmode = 5; %Number of modal coordinates acquiring
ETA = SCORE(:,1:nmode)'; % Get first nmode modal coord.
U_r = COEFF(:,1:nmode); % Get first nmode modal coord.
% ------------------------------------------------------------------------
%LOW-FPS - DOWNSAMPLED(DS)
[COEFFDS, SCOREDS, ~, ~, explainedDS] = pca(Phase_2D_DS,'Algorithm','svd');
%nmode = 5; %Number of modal coordinates acquiring
ETA_DS = SCOREDS(:,1:nmode)'; % Get first nmode modal coord.
U_r_DS = COEFFDS(:,1:nmode); % Get first nmode modal coord.
fprintf('Done\n');

fprintf('Running complexity pursuit algorithm...\n');
% ------------------------------------------------------------------------
%HI-FPS
[A W] = CP(ETA);  % A = Mixing matrix for finding mode shapes (cols=PC)
                    % W = Unmixing matrix for sperating out modes (cols=PC)
% ------------------------------------------------------------------------
%LOW-FPS - DOWNSAMPLED(DS)
[A_DS W_DS] = CP(ETA_DS);
fprintf('Done\n');

fprintf('Building mode shape and modal coordinate matrices...\n');
% ------------------------------------------------------------------------
%HI-FPS
q_t = ETA'*W;     % Q_t is a matrix of modal coordinates
Mode_shape = U_r*A; % Mode_shape is a matrix of all mode shapes
                    % Each column is a full-field mode shape
% ------------------------------------------------------------------------
%LOW-FPS - DOWNSAMPLED(DS)
q_t_DS = ETA_DS'*W_DS;
Mode_shape_DS = U_r_DS*A_DS;
fprintf('Done\n');  

%% Extract frequencies and damping ratios from identified modes
fprintf('Extract frequencies and damping ratios from identified modes...\n');
% HI
dt = 1/v.FrameRate;
t = 0:dt:v.Duration-1/v.FrameRate;
% LO
dt_DS = 1/(v.FrameRate/DS_factor);
t_DS = 0:dt_DS:(v.Duration/DS_factor)-1/(v.FrameRate/DS_factor);

% The time sequences are now the columns of q_t
% HI
N = size(q_t,1);
% LO
N_DS = size(q_t_DS,1);

% Put modal coordinate matrix into fourier domain to extract freq.
% and access only first portion (forget the negatives frequencies)
% HI
Qfft = fft(q_t,[],1); % FFT across rows of q_t
Qfft = Qfft([1:(N/2+1)],:,:);
% LO
Qfft_DS = fft(q_t_DS,[],1);
Qfft_DS = Qfft_DS([1:(N_DS/2+1)],:,:);


% Do same for principal components
% HI
ETA_fft = fft(ETA,[],2); % FFT across cols of SCORE
ETA_fft = ETA_fft(:,[1:(N/2+1)],:);
% LO
ETA_DS_fft = fft(ETA_DS,[],2);
ETA_DS_fft = ETA_DS_fft(:,[1:(N_DS/2+1)],:);

% Create Frequency vector in rad/s
% HI
w0 = 2*pi/(N*dt);
ws = [0:w0:w0*N/2].';
% LO
w0_DS = 2*pi/(N_DS*dt_DS);
ws_DS = [0:w0_DS:w0_DS*N_DS/2].';

% Filter out frequencies below 1st mode (7.2 Hz)
thsh = 0; %min Hz frequency
% HI
Qfft(ws/2/pi<=thsh,:)=0;
% LO
Qfft_DS(ws_DS/2/pi<=thsh,:)=0;

% HI
for ll = 1:size(Qfft,2)
    %Find frequency peaks at each mode
    [~,index] = max(abs(Qfft(:,ll)));
    ws(index)/2/pi;
    %Normalize frequency span
    Q_norm = abs(Qfft(:,ll)/Qfft(index,ll));
    Modal_freq_rad(ll) = ws(index); % Frequencies of each mode
    
%     % Data to left of peak:
%     QL = Q_norm(index:-1:1);
%     wl = ws(index:-1:1);
%     id1=find(QL<=1/sqrt(2),1);
%     w1 = interp1(QL(1:id1),wl(1:id1),1/sqrt(2));
% 
%     % Data to right of peak:
%     QR = Q_norm(index:end);
%     wr = ws(index:end);
%     id2 = find(QR<=1/sqrt(2),1);
%     w2 = interp1(QR(1:id2),wr(1:id2),1/sqrt(2));
% 
%     zt(ll) = (w2-w1)/(2*Modal_freq_rad(ll)); %Damping Ratios
end
Modal_freq_hz = Modal_freq_rad/2/pi;

% LO
for ll = 1:size(Qfft_DS,2)
    %Find frequency peaks at each mode
    [~,index] = max(abs(Qfft_DS(:,ll)));
    ws_DS(index)/2/pi;
    %Normalize frequency span
    Q_norm_DS = abs(Qfft_DS(:,ll)/Qfft_DS(index,ll));
    Modal_freq_rad_DS(ll) = ws_DS(index); % Frequencies of each mode
    
%     % Data to left of peak:
%     QL_DS = Q_norm_DS(index:-1:1);
%     wl_DS = ws_DS(index:-1:1);
%     id1=find(QL<=1/sqrt(2),1);
%     w1_DS = interp1(QL_DS(1:id1),wl_DS(1:id1),1/sqrt(2));
% 
%     % Data to right of peak:
%     QR_DS = Q_norm_DS(index:end);
%     wr_DS = ws_DS(index:end);
%     id2 = find(QR_DS<=1/sqrt(2),1);
%     w2_DS = interp1(QR_DS(1:id2),wr_DS(1:id2),1/sqrt(2));
% 
%     zt_DS(ll) = (w2_DS-w1_DS)/(2*Modal_freq_rad_DS(ll)); %Damping Ratios
end
Modal_freq_hz_DS = Modal_freq_rad_DS/2/pi;

fprintf('Done\n');

%% Extract frequencies from mic data

figure;
%freq = mic.Slot5_ai2.Fs*[0:(mic.Slot5_ai2.L-1)]/(mic.Slot5_ai2.L);
%fft_mic = fftshift(fft(mic.Slot5_ai2.data));
dt = mic.t(1,2)-mic.t(1,1);
Fs = 1/dt;
L = size(mic.t,2);
freq = Fs*[0:(L-1)]/L;
whole_fft_mic = fftshift(fft(mic.x));
fourier_mic = abs(whole_fft_mic(length(whole_fft_mic)/2:end-1));
plot(freq(1:length(freq)/2),fourier_mic); 
title('Frequencies of Mic Data'); xlabel('Frequency (Hz)');

numFreq = nmode;
[max_mic_magn, ind] = maxk(fourier_mic,numFreq);    %Find numFreq number of max 
                                                %values in mic data
for ii=1:numFreq
    fprintf('Frequency');
    freq(ind(1,ii));
    fprintf('Amplitude');
    max_mic_magn(1,ii);
end

%% Perform anti-aliasing to Down-Sampled camera data
FPS_DS_video = round(size(vid_downsample,4)/(L*dt));


%% Extract Mode Shapes
fprintf('Extract mode shape information...\n');
% ------------------------------------------------------------------------
%HI-FPS
ModeShapes = zeros(v.Height,v.Width,nmode,'single');
for uu = 1:nmode
   for oo = 1:size(coord_list,1)
       ii = coord_list(oo,1);
       jj = coord_list(oo,2);
       
       %Fill in the non-noise values
       ModeShapes(ii,jj,uu) = Mode_shape(oo,uu);
   end
end
% ------------------------------------------------------------------------
%LOW-FPS - DOWNSAMPLED(DS)
ModeShapesDS = zeros(v.Height,v.Width,nmode,'single');
for uu = 1:nmode
   for oo = 1:size(coord_listDS,1)
       ii = coord_listDS(oo,1);
       jj = coord_listDS(oo,2);
       
       %Fill in the non-noise values
       ModeShapesDS(ii,jj,uu) = Mode_shape_DS(oo,uu);
   end
end
fprintf('Done\n');

%% Graphical Outputs 

fprintf('Plotting spectrogram of frequencies in video\n');
%%%%%%%%%%%%%%%%%%
%%% Katie Code %%%
%%%%%%%%%%%%%%%%%%
fprintf('Done\n');

fprintf('Plotting CP/PCA modal coordinate graphs, and Mode Shape...\n');
% ------------------------------------------------------------------------
%HI-FPS
figure('Hi-FPS Video');
for dd = 1:nmode
    subplot(nmode,2,2*dd-1);
    plot(t,ETA(dd,:));
    if dd==1
        title('Modes','FontSize',20);
    elseif dd==2||(dd==nmode&&nmode==2)
        ylabel('PCA Coordinate','FontSize',16);
    elseif dd==nmode
        xlabel('Time (s)','FontSize',16);
    end
    subplot(nmode,2,2*dd)
    plot(freq,abs(ETA_fft(dd,:)).^2);
    if dd==1
        title('Power Spectral Density','FontSize',20);
    elseif dd==2||(dd==nmode&&nmode==2)
        ylabel('Amplitude (squared to minimize noise)','FontSize',16);
    elseif dd==nmode
        xlabel('Frequency (Hz)','FontSize',16);
    end
end
figure('Hi-FPS Video');
for ee = 1:nmode
    ModeShapes(:,:,ee)=ModeShapes(:,:,ee)/max(max(abs(ModeShapes(:,:,ee))));
    figure;
    ampl_factor = 100;
    surf(ampl_factor*ModeShapes(:,:,ee)); colormap default; shading interp; 
    axis equal off;
    if ee==1
        title('Mode Shape 3');
    elseif ee==2
        title('Mode Shape 2');
    else
        title('Mode Shape 1');
    end
end
% ------------------------------------------------------------------------
%LOW-FPS - DOWNSAMPLED(DS)
figure('Low-FPS Video');
for dd = 1:nmode
    subplot(nmode,2,2*dd-1);
    plot(t,ETA(dd,:));
    if dd==1
        title('Modes','FontSize',20);
    elseif dd==2||(dd==nmode&&nmode==2)
        ylabel('PCA Coordinate','FontSize',16);
    elseif dd==nmode
        xlabel('Time (s)','FontSize',16);
    end
    subplot(nmode,2,2*dd)
    plot(freq,abs(ETA_fft(dd,:)).^2);
    if dd==1
        title('Power Spectral Density','FontSize',20);
    elseif dd==2||(dd==nmode&&nmode==2)
        ylabel('Amplitude (squared to minimize noise)','FontSize',16);
    elseif dd==nmode
        xlabel('Frequency (Hz)','FontSize',16);
    end
end
figure('Low-FPS Video');
for ee = 1:nmode
    ModeShapes(:,:,ee)=ModeShapes(:,:,ee)/max(max(abs(ModeShapes(:,:,ee))));
    figure;
    ampl_factor = 100;
    surf(ampl_factor*ModeShapes(:,:,ee)); colormap default; shading interp; 
    axis equal off;
    if ee==1
        title('Mode Shape 3');
    elseif ee==2
        title('Mode Shape 2');
    else
        title('Mode Shape 1');
    end
end
fprintf('Done\n');
