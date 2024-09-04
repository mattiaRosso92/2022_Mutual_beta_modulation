function [filtdat] = filter_plateau(signal,srate,frange,trans_width, order , showplot)
%  filterFGx   Narrow-band filter via plateau-shaped FIR filter
%  
%    Reference: https://www.udemy.com/course/solved-challenges-ants/learn/lecture/13322022#overview
% 
%    INPUTS
%       data : 1 X time or chans X time
%      srate : sampling rate in Hz
%     frange : vector with boundary frequencies (e.g. [8 12]) 
% trans_width : transation width (i.e., slope in frequency domain); in
%              decimals (e.g., .2 = 20%)
%   showplot : set to true to show the frequency-domain filter shape
% 
%    OUTPUTS
%    filtdat : filtered data
% 
% 
% mikexcohen@gmail.com

%% input check

if size(signal,1)>size(signal,2)
%     help filter_plateau
%     error('Check data size')
end

if nargin<5
    help filter_plateau
    error('Not enough inputs')
end

if frange(end)<=frange(1)
    error('Insert valid frequency range')
end

if nargin<6
    showplot=false;
end

%% Compute filter

% Settings
nyquist = srate/2;
filt_order   = round(order*srate/frange(1)); %length of the kernel in time domain (i.e., N of points); 
% the longer the filter, the higher frex res; it also slows down the
% analyses, and if exaggerated it would give weird effects in frex domain;
% check that the frequency response is nicely flat at 1 in the desired
% range
% I'd say take advantage of the long recordings to use very high order
% filter and boost frequency resolution

% 6 points of the filter response (as fraction of nyquist)
frex_vector = [0 (1-trans_width)*frange(1) frange (1+trans_width)*frange(2) nyquist] / nyquist; 
% ideal response (same number of points as ffrequencies)
ideal_response = [0 0 1 1 0 0];

%get filter weigths
filter_weights = firls(filt_order,frex_vector,ideal_response);

% Create FIR filter

% Filter Kernel
%filtkern = fir1(order,frange/nyquist); %general function, regardless of unit and domain: second argument in fraction of the nyquist
% The slope is given by default; check-out the more general firls()
% function for modifying it

% % Compute the power spectrum of the filter (i.e., filter in frex domain)
% filtpow = abs(fft(filtkern)).^2;
% % Compute frequencies vector and remove negative frex
% hz = linspace(0,srate/2,floor(length(filtkern)/2)+1);
% filtpow = filtpow(1:length(hz));

%% Filter the input signal
% 
% % Initialize output
filtdat = zeros(size(signal));
% Filter channel-by-channel 
for chani = 1:size(signal,1) 
    
    filtdat(chani,:) = filtfilt(filter_weights,1,signal(chani,:)); %filtfilt() corrects phase shifts

end

%(maybe reshape and filter long)
% nchans = size(data,1);
% data = reshape(data,1,[]);
% filtdat  = filtfilt(filtkern,1,data);
% filtdat = reshape(filtdat,nchans,[]);


%% inspect the FIR filter (turned off by default)
% 
% if showplot
%     figure(10001+showplot),clf
%     sgtitle('Required filter')
%     % Filter in time domain
%     subplot(121)
%     plot(filtkern,'k','linew',2)
%     xlabel('Timepoints')
%     title('Filter Kernel (fir1)')
%     axis square
%     % Filter in frequency domain
%     subplot(122) , hold on
%     plot(hz,filtpow,'ks-','linew',2,'MarkerFaceColor','w') % actual filter (slopes are given by fir1())
%     plot([0 frange(1) frange frange(end) nyquist],[0 0 1 1 0 0],'ro-','linew',1.5)% ideal filter
%     set(gca,'xlim',[frange(1)-3 frange(end)+3])
%     xlabel('Frequency (Hz)'), ylabel('Filter gain')
%     title('Frequency response of filter')
%     legend('Actual','Ideal')
%     axis square
% 
% end

if showplot
    figure(10001+showplot),clf
%     sgtitle('Requested filter')
    % Filter in frequency domain
%     subplot(211)
    plot(frex_vector*nyquist,ideal_response,'k-o','MarkerFace','m', 'LineWidth', 1.3)
    % TODO: PLOT ACTUAL RESPONSE OVER REQUESTED; CHECK THEY MATCH
    set(gca,'ylim',[-.1 1.1],'xlim',[0 frange(end)+7])
    xlabel('Frequencies (Hz)')
    ylabel('Gain')
    title('Frequency response of filter')
%     % Filter in time domain
%     subplot(212) , hold on
%     plot( (0:filt_order) * (1000/srate) , filter_weights , 'k' , 'linew' , 1.3)
%     xlabel('Time (ms)'), ylabel('Amplitude')
%     title('Filter kernel')
%     
   
end

%% done.
