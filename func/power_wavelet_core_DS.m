function [P1,P2,freqs2use,t]= power_wavelet_core_DS(s_1,s_2)
%% Downsampling to 200 Hz
fs=200;
t=decimate([1:2100],5);
s_1_= []; s_2_ = [];
for trii=1:size(s_1,1)
    s_1_=[s_1_ [ zeros(1,500) (s_1(trii,:)) zeros(1,500)]];
    s_2_=[s_2_ [ zeros(1,500) (s_2(trii,:)) zeros(1,500)]];
     
%      s_1_=[s_1_ (s_1(trii,:))];
%      s_2_=[s_2_ (s_2(trii,:))];
end
s_1_=decimate(s_1_,5,'fir');
s_2_=decimate(s_2_,5,'fir');
freqs2use  =logspace(log10(1),log10(100),60); % 4-100 Hz in 30 steps
%  freqs2use  =1:100; % 4-100 Hz in 30 steps

timewindow = linspace(1,3,length(freqs2use)); % number of cycles on either end of the center point (1 means a total of 2 cycles))
pnts=size(s_1,2)/5+200;trials=size(s_1,1);
% wavelet and FFT parameters
time          = -1:1/fs:1;
half_wavelet  = (length(time)-1)/2;
num_cycles    = logspace(log10(4),log10(8),length(freqs2use));
num_cycles    = 7*ones(1,length(freqs2use));

n_wavelet     = length(time);
n_data        = pnts*trials;
n_convolution = n_wavelet+n_data-1;

% data FFTs
data_fft1 = fft(reshape(s_1_,1,n_data),n_convolution);
data_fft2 = fft(reshape(s_2_,1,n_data),n_convolution);

%% core 
for fi=1:length(freqs2use)
    
    % create wavelet and take FFT
    s = num_cycles(fi)/(2*pi*freqs2use(fi));
    wavelet_fft = fft( exp(2*1i*pi*freqs2use(fi).*time) .* exp(-time.^2./(2*(s^2))) ,n_convolution);
    
    % phase angles from channel 1 via convolution
    convolution_result_fft = ifft(wavelet_fft.*data_fft1,n_convolution);
    convolution_result_fft = convolution_result_fft(half_wavelet+1:end-half_wavelet);
    power_sig1 = (reshape(convolution_result_fft,pnts,trials));
    
    % phase angles from channel 2 via convolution
    convolution_result_fft = ifft(wavelet_fft.*data_fft2,n_convolution);
    convolution_result_fft = convolution_result_fft(half_wavelet+1:end-half_wavelet);
    power_sig2 = (reshape(convolution_result_fft,pnts,trials));
    
    % phase angle differences
    %  phase_diffs = phase_sig1-phase_sig2;
    P1(fi,:)=nanmean(abs(power_sig1).^2,2);
    P2(fi,:)=nanmean(abs(power_sig2).^2,2);
   
    %             end
end % end frequency loop

P1=P1(:,101:520);
P2=P2(:,101:520);
end
