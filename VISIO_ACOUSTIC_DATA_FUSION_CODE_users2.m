% ------------------------------------------------------------------------
% VISIO-ACOUSTIC DATA FUSION FOR STRUCTURAL HEALTH MONITORING APPLICATIONS 
% ------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%% USER EDITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%       Low-FPS (subNyquist) video data
%       Hi-Fs microphone data (from 1 microphone)
% OUTPUT:
%       Spectrogram of frequencies in video
%       Location of damaged object in scene emitting hi-frequency sound
%       Image of scene with outlines around each area per frequency
% ------------------------------------------------------------------------
clc
clear all
close all

%% Load data
% Microphone file
[mic_name,mic_folder]=uigetfile(fullfile(pwd,'folder1','*m4a'),...
    'Select a MICROPHONE file');
if mic_name==0 
    disp('No file selected, exiting...'); 
    return; 
end
[mic Sr] = audioread(mic_name);
L_mic=length(mic);
fprintf('Microphone data loaded\n');
%%
% Video file
[vid_name,vid_folder]=uigetfile(fullfile(pwd,'folder1','*.mov;*.avi;*mp4'),...
    'Select a VIDEO file');
if vid_name==0 
    disp('No file selected, exiting...'); 
    return; 
end

v = VideoReader(vid_name);
% frame=readFrame(v);
% figure;
% subplot(2,1,1);
% imshow(frame);
% subplot(2,1,2);
% imshow(frame(:,v.Width/4+50:v.Width/2));

% Set frame rate
answer=inputdlg({'Enter the frame rate of your video (in frame per second):'},...
    'User Input',1,{num2str(v.FrameRate),'center'});
framerate=str2double(answer{1});

h=waitbar(0,'Loading Video');
if framerate>100
    v.Current=v.Duration/1.5;
    vid=zeros(v.Height,v.Width/2-300,3,floor(v.Duration*v.FrameRate));
else
    v.Current=0;
    vid=zeros(v.Height,v.Width/2-300,3,floor(v.Duration*v.FrameRate));
end

nFrames = 1;
while hasFrame(v)
    frame = readFrame(v);
    %cropframe=frame(:,v.Width/4-100:v.Width/2-1,:); %60fps v3 - 512Hz
    %cropframe=frame(:,v.Width/2+100:v.Width-1,:); %60fps v3 - 128Hz
    %cropframe=frame(:,v.Width/4:v.Width/2-1,:); %120fps v3 - 512Hz
    cropframe=frame(:,v.Width/2+200:v.Width-100-1,:); %120fps v3 - 128Hz
    %cropframe=frame(:,v.Width/4+50:v.Width/2-1,:); %240fps v3 - 512Hz
    %cropframe=frame(:,v.Width/2+100:v.Width-200-1,:); %240 v3 - 128Hz
    vid(:,:,:,nFrames) = cropframe; 
    waitbar(v.CurrentTime/v.Duration,h)
    nFrames = nFrames+1;
    clear frame;
end
nFrames = nFrames-1;
close(h); 
fprintf('Video Data Loaded\n');

%% Filter video data, build 3D & 2D Phase and Amplitude Matrices
fprintf('Filtering video and building Phase & Amplitude Matrices\n');

ht = maxSCFpyrHt(zeros(v.Height,v.Width/2-300));
pyrType = 'octave';
switch pyrType
    case 'octave'
        filters = getFilters([v.Height v.Width/2-300], 2.^[0:-1:-ht], 4);
        repString = 'octave';
        fprintf('Using octave bandwidth pyramid\n');   
    otherwise 
        error('Invalid Filter Types');
end
vidFFT = zeros(v.Height,v.Width/2-300,nFrames,'single');

[croppedFilters, filtIDX] = getFilterIDX(filters);

buildLevel = @(im_dft, k) ifft2(ifftshift(croppedFilters{k}.* ...
        im_dft(filtIDX{k,1}, filtIDX{k,2})));

fprintf('Moving video to Fourier domain\n');
h=waitbar(0,'Moving video to Fourier domain');
for k = 1:nFrames
    originalFrame = rgb2ntsc(im2single(vid(:,:,:,k)));
    tVid = imresize(originalFrame(:,:,1), [v.Height v.Width/2-300]);
    vidFFT(:,:,k) = single(fftshift(fft2(tVid)));
    waitbar(k/nFrames,h);
end
close(h);

level = 2; %1=Hi; 2=vert. filt(-); 3=45 cw(\); 4=horiz. filt(l); 5=135 cw(/)...

% Obtain reference phase/frame data
pyrRef = buildLevel(vidFFT(:,:,1), level);
%pyrRefPhaseOrig = pyrRef./abs(pyrRef);
phaseRef = angle(pyrRef);

fprintf('Building 3D Amplitude & Phase matrices\n');
Phase = zeros(size(pyrRef,1), size(pyrRef,2) ,nFrames,'single');
Amp = zeros(size(pyrRef,1), size(pyrRef,2) ,nFrames);
for frameIDX = 1:nFrames
    filterResponse = buildLevel(vidFFT(:,:,frameIDX), level);
    pyrCurrent = angle(filterResponse);
    Amp(:,:,frameIDX) = abs(filterResponse);
    Phase(:,:,frameIDX) = single(mod(pi+pyrCurrent-phaseRef,2*pi)-pi);                          
end
fprintf('Done\n\n');

%% Anti-Aliasing Idea
fprintf('Perform Anti-Aliasing Algorithm to DownSampled video data\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% ANTI-ALIASING STEPS %%%%%%%%%%%%%%%%%%%%%%%%%
% 1. DOWN SAMPLE MICROPHONE DATA
%       1.a)    Filter out frequencies below 100 Hz on orig mic data
%       1.b)    Down Sample Mic data
%       1.c)    Convert Mic data into fourier domain
% 2. OVERLAY AND MATCH FFT PLOTS
%       2.a)    Use findpeaks function to find peaks in mic data
% 3. Make DS Alias calculation & Filter DS video accordingly
%       3.a)    Create alias finding calculation on found peak frequencies
%       3.b)    Plot peak frequencies found in mic data and aliased
%               frequencies in lollipop plotting format
%       3.c)    Filter DS video data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 1.a: Filter out frequencies below 100 Hz
micfilt=highpass(mic,100,Sr); 

%% STEP 1.b: DOWN SAMPLE MICROPHONE DATA
fprintf('Down Sampling Microphone Data\n');
% Subsample mic data
micDS = zeros(floor(L_mic/Sr/framerate),1);
Length_DSmic = 1; %number of points in subnyquist mic data
for ii=1:floor(Sr/framerate):L_mic
    micDS(Length_DSmic,1) = micfilt(ii,1);
    Length_DSmic=Length_DSmic+1;
end
Length_DSmic=Length_DSmic-1;

%% STEP 1.c: Convert Mic data into fourier domain
fprintf('Convert original & down sampled mic data into fourier domain\n');
%Transform mic to fourier domain
micfreq = Sr*[0:(L_mic-1)]/(L_mic);
micFFT = fft(micfilt);
micfreq_pos = micfreq(1:length(micfreq)/2);
micFFT_pos = abs(micFFT(1:length(micFFT)/2));
micDSfreq = framerate*[0:(Length_DSmic-1)]/Length_DSmic;
micDSFFT = fft(micDS);
micDSfreq_pos = micDSfreq(1:length(micDSfreq)/2);
micDSFFT_pos = abs(micDSFFT(1:length(micDSFFT)/2));

%% STEP 2: OVERLAY AND MATCH FFT PLOTS
%% STEP 2.a: Use findpeaks function to find peaks in mic data
%% Mic data
fprintf('Finding peaks in orginal mic data\n');
% Original Mic data
[PKSmic,location_mic_Hz]=findpeaks(micFFT_pos,micfreq_pos,'MinPeakHeight',.2,...
    'MinPeakDistance',20,'SortStr','descend','NPeaks',3);
PKSmic;
location_mic_Hz;

%% STEP 3: Make DS Alias calculation & Filter DS video accordingly
%% STEP 3.a:    Create alias finding calculation on found peak frequencies
alias_freqs=zeros(1,length(location_mic_Hz));
for ii=1:length(location_mic_Hz)
    currFreq=location_mic_Hz(ii); %Frequency value in Hz
    currDec=currFreq/(framerate/2); %Ratio of the Frequency value to 1/2 sampling rate
    currInt=floor(currFreq/(framerate/2)); %Same ratio as above, but decimal is truncated
    if mod(currInt,2)==0 %Even
        remaining=currDec-currInt;  %remaining decimal portion represents 
                                    %fraction of 1/2 the sampling rate
        alias_freqs(ii)=remaining*(framerate/2);
    else %Odd
        remaining=currDec-currInt;  %remaining decimal portion represents 
                                    %1-fraction of 1/2 the sampling rate
        alias_freqs(ii)=(1-remaining)*(framerate/2);
    end
end
lolivec=zeros(1,length(micfreq_pos));
for ii=1:length(location_mic_Hz)
    [val idx]=min(abs(micfreq_pos-alias_freqs(ii)));
    lolivec(idx)=PKSmic(ii);
end
fprintf('Corresponding aliased frequencies calculated\n');

%% STEP 3.b:    Plot found peaks against aliased peaks (lollipop format)
figure;
subplot(2,1,1);
findpeaks(micFFT_pos,micfreq_pos,'MinPeakHeight',.2,...
    'MinPeakDistance',20,'SortStr','descend','NPeaks',3,'Annotate','peaks');
title(['Peaks found in original mic data(>20 Hz apart) @ ',num2str(Sr),' Sampling Rate'],'FontSize',20);
xlabel('Frequency (Hz)','FontSize',20);  ylabel('|[a.u.](f)|','FontSize',20);
xticks([1:5:length(micfreq_pos)]*100); xtickangle(45);
text(location_mic_Hz+.02,PKSmic,num2str((1:numel(PKSmic))'))
set(gca,'FontSize',16);
legend('Original Microphone Data','Main 3 peaks found');

subplot(2,1,2);
stem(micfreq_pos,lolivec); title(['Found peaks aliased @ ',num2str(framerate),' Sampling Rate'],'FontSize',20);
xlabel('Frequency (Hz)','FontSize',20);  ylabel('|[a.u.](f)|','FontSize',20);
xticks([1:5:length(micfreq_pos)]*100); xtickangle(45);
text(alias_freqs-.5,PKSmic+.1,num2str((1:numel(PKSmic))'))
set(gca,'FontSize',16); 
hold on;
xline(framerate/2,'--r');
legend('Aliased folded frequencies',['Subnyquist Frequency (',num2str(framerate/2),')']);

figure;
findpeaks(micFFT_pos,micfreq_pos,'MinPeakHeight',.2,...
    'MinPeakDistance',20,'SortStr','descend','NPeaks',3,'Annotate','peaks');
title(['Aliased Frequency Peaks overlayed on top of Original Microphone Data Peaks with Subnyquist Frequency @ ',...
    num2str(framerate),' Sampling Rate'],'FontSize',20);
xlabel('Frequency (Hz)','FontSize',20);  ylabel('|[a.u.](f)|','FontSize',20);
xticks([1:5:length(micfreq_pos)]*100); xtickangle(45);
text(location_mic_Hz+.02,PKSmic,num2str((1:numel(PKSmic))'))
set(gca,'FontSize',16); hold on;

stem(micfreq_pos,lolivec); 
text(alias_freqs-.5,PKSmic+.01,num2str((1:numel(PKSmic))'))
hold on;
xline(framerate/2,'--g');
legend('Original Microphone Data','Main 3 peaks found','Aliased folded frequencies',...
    ['Subnyquist Frequency (',num2str(framerate/2),')']);

%% STEP 3.c:    Bandpass Filter video data based on aliased mic data peaks
fprintf('Bandpass filtering video data based on aliased frequencies\n');
% Using temporal filter function
temporalFilter = @FIRWindowBP;
delta_filt4D=zeros(size(Phase,1),size(Phase,2),size(Phase,3),length(location_mic_Hz),'single');
fs=framerate;
h=waitbar(0,'Filtering video');
for ii=1:length(location_mic_Hz)
    %Bandpass filter around alias frequency
    if (alias_freqs(ii)-5)>=0
        fl=alias_freqs(ii)-2;
    else
        fl=alias_freqs(ii);
    end
    if (alias_freqs(ii)+5)<(fs/2)
        fh=alias_freqs(ii)+2;
    else
        fh=alias_freqs(ii);
    end
    delta_filt4D(:,:,:,ii) = temporalFilter(Phase, fl/fs,fh/fs);
    waitbar(ii/length(location_mic_Hz),h);
end
close(h);

fprintf('Formatting filtered video to be used for Segmentation\n');
magPhase=200;
sigma=5;
for kk=1:length(location_mic_Hz)
    for ii=1:nFrames
        delta_new(:,:) = delta_filt4D(:,:,ii,kk);
        % Amplitude Weighted Blur        
        if (sigma~= 0)
            delta_smooth = AmplitudeWeightedBlur(delta_new, Amp(:,:,ii)+eps, sigma);        
        end
        % Increase phase variation
        delta_smooth = magPhase *delta_smooth;  
        delta_filt4D(:,:,ii,kk) = magPhase*delta_smooth;
    end
end
%delta_filt3D=[pixel_count,frame_num,filtered_freqs];
delta_filt3D=reshape(delta_filt4D,[size(Phase,1)*size(Phase,2),size(Phase,3),length(location_mic_Hz)]);
fprintf('Done\n\n');
Amp1=Amp(:,:,1);
vid1=vid(:,:,:,1); 
%% PERFORM SEGMENTATION ON FILTERED DATA SET
fprintf('Performing Segmentation on Video\n');
close all
clc
%inFile='tuningforks_128_512_2048Hz_93in_5000sr.mov';
%Ptime_xy=zeros(size(vid,1),size(vid,2),length(location_mic_Hz));
rgb=zeros(3480,size(delta_filt4D,2),3,length(location_mic_Hz));
yellowPixels=zeros(1,length(location_mic_Hz));
for ii=1:length(location_mic_Hz)
power=zeros(size(delta_filt3D,1),size(delta_filt3D,2),length(location_mic_Hz));
sd=zeros(size(delta_filt3D,1),1,length(location_mic_Hz));
 
for jj=1:size(delta_filt3D,1)
%Finding power of each pixel across all frequencies 
power(jj,:,ii)=abs(fftshift(fft2(delta_filt3D(jj,:,ii)),2)).^2;
sd(jj,1,ii)=sum(power(jj,:,ii));
end
%Reshapes sd to replicate matrix of video (i.e. x vs y coordinates) 
Ptime_xy=reshape(sd,size(delta_filt4D,1),size(delta_filt4D,2),length(location_mic_Hz));

% fig=imagesc(Ptime_xy(:,:,ii));

set(gca,'YTickLabel',[])
set(gca,'XTickLabel',[])
% 
%  saveas(fig, sprintf('128Hzimagesc%i.png',ii));
% 
% BB(:,:,:,ii)=imread(sprintf('128Hzimagesc%i.png',ii));
%  
% lab_he = rgb2lab(BB(:,:,:,ii));
%  ab = lab_he(:,:,1);
% ab = im2single(ab);
nColors = 3;


pixel_labels = imsegkmeans(single(Ptime_xy(:,:,ii)),nColors,'NumAttempts',3);
mask3 = pixel_labels==3;
cluster1 = Ptime_xy(:,:,ii) .* mask3;
imshow(cluster1)


%ww=imagesc(cluster1);
%hold on

%I=imbinarize(rgb2gray(cluster1));

er=edge(cluster1);
[row,col]=find(er);

row_new=zeros(size(row,1),1);
col_new=zeros(size(row,1),1);
row_new2=zeros(size(row,1),1);
col_new2=zeros(size(row,1),1);
for j=1:length(row)
    row_new(j)=row(j);
      col_new(j)=col(j);
       row_new2(j)=row(j);
      col_new2(j)=col(j);
end 
row_new=cat(1,row_new,row,row_new2);
col_new=cat(1,col_new,col,col_new2);

origframe=imagesc(vid(:,:,:,1));
set(gca,'YTickLabel',[])
set(gca,'XTickLabel',[])


% saveas(origframe,'firstframe3tuning_forks.png');
%  A=imread('firstframe3tuning_forks.png');
 A=origframe.CData;
% A=A.*0.25;
A=zeros(size(origframe.CData,1),size(origframe.CData,2),size(origframe.CData,3));
r=[250,0,0,250,0,125,125];
g=[0,250,0,0,250,125,50];
b=[0,0,250,250,125,0,50];
for jj=1:length(row_new)
    A(row_new(jj),col_new(jj),:)=cat(3,[r(ii)],[g(ii)],[b(ii)]);
end
xx=imagesc(A);
AA(:,:,:,ii)=xx.CData;
yy=imagesc(Ptime_xy(:,:,ii));
BB(:,:,:,ii)=yy.CData;
cm = colormap(yy.Parent); % get axes colormap
n = size(cm,1); % number of colors in colormap
c = linspace(yy.Parent.CLim(1),yy.Parent.CLim(2),n); % intensity range
ind = reshape(interp1(c,1:n,Ptime_xy(:,:,ii),'nearest'),size(Ptime_xy(:,:,ii))); % indexed image
rgb(:,:,:,ii) = ind2rgb(ind,cm); % rgb image
figure;
imshow(rgb(:,:,:,ii),'InitialMagnification','fit');

yellowPixels(ii)= size((find(rgb(:,:,1,ii) >0.75 & rgb(:,:,2,ii) > .75 & rgb(:,:,3,ii) < 0.25)),1);

% numYellowPixels(ii) = sum(yellowPixels(ii),'all');

close all
end
freq_edges=zeros(size(AA,1),size(AA,2),size(AA,3),size(AA,4));
freq_edges=AA(:,:,:,1);
for ii=1:length(location_mic_Hz)
    freq_edges=freq_edges+AA(:,:,:,ii);
end
freq_amps=zeros(size(rgb,1),size(rgb,2),size(rgb,3),size(rgb,4));
index=find(yellowPixels>mean(yellowPixels));
for ii=index
    freq_amps=rgb(:,:,:,index(1));
    freq_amps=freq_amps+rgb(:,:,:,ii);
%    imshow(freq_amps);
end
[rowblue,colblue]=find(freq_amps(:,:,1)<.5 & freq_amps(:,:,2) <0.32 & freq_amps(:,:,3) > 1.32);

for jj=1:length(rowblue)
    freq_amps(rowblue(jj),colblue(jj),:)=cat(3,0,0,0);
end
fused_segmented_image=0.7.*freq_amps+freq_edges+.1.*Amp(:,:,1);
%% Show Segmented Image
figure;
imshow(fused_segmented_image); title('Segmented Image','FontSize',20);
text(1, 230,'Segmented','BackgroundColor','w','Color','k','FontSize',15,'FontWeight','Bold','Interpreter','None');
text(1, 250,'Frequencies (Hz) ±5 Hz:','BackgroundColor','w','Color','k','FontSize',15,'FontWeight','Bold','Interpreter','None');
for ii=1:length(location_mic_Hz)
    str = sprintf('%i Hz',round(location_mic_Hz(ii)));
    text(1, 250+ii*20,str,'BackgroundColor','w','Color',[r(ii)/250 g(ii)/250 b(ii)/250],'FontSize',15,'FontWeight','Bold','Interpreter','None');
end
%% User friendly Segmented Image:
f=figure;
p=uipanel(f,'Position',[0 0 .15 .125]);
c=uicontrol(p,'Style','popupmenu','FontSize',20);
c.Value=1;
c.String={'Combined Frequency Segmented Image',sprintf('Only %i Hz',round(location_mic_Hz(1))),...
    sprintf('Only %i Hz',round(location_mic_Hz(2))),sprintf('Only %i Hz',round(location_mic_Hz(3))),...
    'Show Original Image'};
c.Callback={@selection,fused_segmented_image,location_mic_Hz,r,g,b,rgb,AA,Amp};

function selection(src,~,fused_segmented_image,location_mic_Hz,r,g,b,rgb,AA,Amp)
    val=src.Value;
    switch val
        case 1
            imshow(fused_segmented_image); title('Segmented Image','FontSize',20);
%             text(1, 230,'Segmented','BackgroundColor','w','Color','k','FontSize',15,'FontWeight','Bold','Interpreter','None');
%             text(1, 250,'Frequencies (Hz) ±5 Hz:','BackgroundColor','w','Color','k','FontSize',15,'FontWeight','Bold','Interpreter','None');
%             for ii=1:length(location_mic_Hz)
%                 str = sprintf('%i Hz',round(location_mic_Hz(ii)));
%                 text(1, 250+ii*20,str,'BackgroundColor','w','Color',[r(ii)/250 g(ii)/250 b(ii)/250],'FontSize',15,'FontWeight','Bold','Interpreter','None');
%             end
        case 2
            imshow(0.7.*rgb(:,:,:,val-1)+AA(:,:,:,val-1)+.1.*Amp(:,:,1));
%             text(1, 230,'Segmented','BackgroundColor','w','Color','k','FontSize',15,'FontWeight','Bold','Interpreter','None');
%             text(1, 250,'Frequencies (Hz) ±5 Hz:','BackgroundColor','w','Color','k','FontSize',15,'FontWeight','Bold','Interpreter','None');
%             str = sprintf('%i Hz',round(location_mic_Hz(val-1)));
%             text(1, 270,str,'BackgroundColor','w','Color',[r(val-1)/250 g(val-1)/250 b(val-1)/250],'FontSize',15,'FontWeight','Bold','Interpreter','None');
        case 3
            imshow(0.7.*rgb(:,:,:,val-1)+AA(:,:,:,val-1)+.1.*Amp(:,:,1));
%             text(1, 230,'Segmented','BackgroundColor','w','Color','k','FontSize',15,'FontWeight','Bold','Interpreter','None');
%             text(1, 250,'Frequencies (Hz) ±5 Hz:','BackgroundColor','w','Color','k','FontSize',15,'FontWeight','Bold','Interpreter','None');
%             str = sprintf('%i Hz',round(location_mic_Hz(val-1)));
%             text(1, 270,str,'BackgroundColor','w','Color',[r(val-1)/250 g(val-1)/250 b(val-1)/250],'FontSize',15,'FontWeight','Bold','Interpreter','None');
        case 4
            imshow(0.7.*rgb(:,:,:,val-1)+AA(:,:,:,val-1)+.1.*Amp(:,:,1));
%             text(1, 230,'Segmented','BackgroundColor','w','Color','k','FontSize',15,'FontWeight','Bold','Interpreter','None');
%             text(1, 250,'Frequencies (Hz) ±5 Hz:','BackgroundColor','w','Color','k','FontSize',15,'FontWeight','Bold','Interpreter','None');
%             str = sprintf('%i Hz',round(location_mic_Hz(val-1)));
%             text(1, 270,str,'BackgroundColor','w','Color',[r(val-1)/250 g(val-1)/250 b(val-1)/250],'FontSize',15,'FontWeight','Bold','Interpreter','None');
        case 5
            imshow(.5.*Amp(:,:,1)); title('Original Image 1st frame','FontSize',20);
        otherwise
            error('this cannot happen');
    end
end