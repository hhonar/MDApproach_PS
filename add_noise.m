function r = add_noise(data,mode_memd,intensity_noise,N_channel_add)
%add_noise() Add white Gaussian noise channel(s) 
%
% Input
%   1) data, signal input formed by a 2D matrix with n rows and m columns.The number of rows indicate the number of data samples for time series of each channel. And the number of columns show the number of channels.
%   2) mode_memd, options for assited noise£º
%     a) mode_memd = 'memd' : no any noise is added, and return 'data';
%     b) mode_memd = 'na_snr' : According to the option of 'na_snr', the noise channels 'r' are returned by tuning the paparmeter 'intensity_noise';
%     c) mode_memd = 'na_fix':  According to the option of 'na_snr', the noise channels 'r' are returned by tuning the paparmeter 'intensity_noise';
%   3) intensity_noise, tuning parameter for the strengths of assisted noise.In 'na_snr' mode, it is the signal-to-noise ratio in dB (the ratio of the power of the signal and the power of assisted Gaussian White noise)£»in 'na_fix' mode, it is the ampltiude of the standard devition of assisted Gaussian White noise.
%   
%   4) N_channel_add, the number of channels for noise signals 
%
% Output
%   r, data with noise channel(s)
%
% Others
%   1) The lengths of noise channels should be equal to those of orignial signal channels;
%   2) If the number of channel in orignial signal is N, and noise channels are M, then the returned output 'r' from 1st to Nth channels presents the corresponding origial signal channels, the remaining from (N+1)th to (N+M)th channels are the assisted noise channels.

r = data;
N_length = size(data,1);

if strcmp(mode_memd, 'na_snr');
    snr_add = intensity_noise;
    sigPower = sum(abs(r(:)).^2)/ N_length;
    noise_add = sqrt(sigPower/(10^(snr_add/10))).* randn(N_length, N_channel_add);
    r = [r, noise_add]; 
elseif strcmp(mode_memd, 'na_fix');
    std_add = intensity_noise;
    noise_add = std_add.*randn(N_length, N_channel_add);
    r = [r, noise_add];
elseif ~strcmp(mode_memd, 'memd');
    Nmsgid = generatemsgid('invalid mode_memd');
    error(Nmsgid,'mode_memd should be a string in na_snr(default) or na_fix');
end
end
