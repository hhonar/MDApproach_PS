% function NAMEMD applies the "Noise Assisted Multivariate Empirical Mode Decomposition" algorithm to multivariate inputs, 
% which can be applied in EEG and EMG signal processing.
% 
% This function extends the "Multivariate Empirical Mode Decomposition" algorithm" (Rehman,Park,Huang, and Mandic,
% Advances in Adaptive Data Analysis, 2013). 
%
% Syntax:
% 
% a)imfs = namemd(x)
%   returns a cell matrix 'imf(N,1)' containing N multivariate Intrinsic Mode Functions (IMFs). 
%   It only can decompose the signal x by using "Multivariate Empirical Mode Decomposition" algorithm. 
%   No noises are added with the signal x for EMD-based decompositions.  
%     - x, original multichannel signals, are formed by a matrix with n rows and m columns.
%       The number of rows indicate the number of data samples for time series of each channel. 
%       And the number of columns show the number of channels. If m > n, the transpose of x is then done 
%       since the number of rows should be far more than the number of columns in the default setting. 
%       The total number of projections of the signal, ndir, is set to 64 (default).      
%     - imfs in the formate of cell. Each component in cell represents the decomposed IMFs in the corresponding channel. 
%
% b)imfs = namemd(x,ndir)
%   where integer variable ndir (>= 1) specifies the total number of projections of the signal, x,
%     - As a rule of thumb, the minimum value of ndir should be twice of the number of data channels,
%     - for instance, ndir = 6  for a 3-variate signal and num_directions= 16 for an 8-variate signal
%   The default number of directions is chosen to be 64 - to extract meaningful IMFs, the number of directions
%   should be considerably greater than the dimensionality of the signals
% 
% c)imfs = namemd(x,ndir,stp_crit)
%   uses the optional parameter stp_crit to control the sifting process.
%    The available options are
%    - 'stop' which uses the standard stopping criterion specified in [2]
%    - 'fix_h' which uses the modified version of the stopping criteria specified in [3]
%    The default value for the 'stopping criteria' is 'stop'.
%  The settings  num_directions=64 and 'stopping criteria' = 'stop' are defaults.
%     Thus imf = MEMD(X) = MEMD(X,64) = MEMD(X,64,'stop') = MEMD(X,[],'stop'),
%
% imf = namemd(X, num_directions, 'stop', stop_vec)
%   computes the IMFs based on the standard stopping criterion whose parameters are given in the 'stop_vec'
%     - stop_vec has three elements specifying the threshold and tolerance values used.
%     - the default value for the stopping vector is   step_vec = [0.075 0.75 0.075].
%     - the option 'stop_vec' is only valid if the parameter 'stopping criteria' is set to 'stop'.
%
% imf = namemd(X, num_directions, 'fix_h', n_iter)
%   computes the IMFs with n_iter (integer variable) specifying the number of consecutive iterations when
%   the number of extrema and the number of zero crossings differ at most by one.
%     - the default value for the parameter n_iter is set to n_iter = 2.
%     - the option n_iter is only valid if the parameter  'stopping criteria' = 'fix_h'
%
% d)imfs = namemd(x, ndir, stp_crit, stp_vec, mode, intensity_noise, n_channel_na)
%   provides the options of specifying the setting of the Gaussian White noise for the noise-assisted analysis of the multichannel signals. 
%   The options can be set by modifying the variable "mode" by means of: 
%     - 'memd' , standard Multivariate Empirical Mode Decomposition (MEMD) mode without noise assistances.It is the default setting.
%     - 'na_snr',the amplitude of the added noise is determined according to the signal-to-noise ratio 
%      (computed by 10lg(PS/PN) where PS and PN indicate the average power of the original signal and added noise). 
%     - 'na_fix', the amplitude of the added noise is determined according to the standard deviation of the Gaussian White noise. 
% 
% imfs = namemd(x, ndir, stp_crit, stp_vec, 'memd')
%   uses the standard MEMD algorithm to decompose the signal x, and obtains the decomposed IMFs. 
%
% imfs = namemd(x, ndir, stp_crit, stp_vec, 'na_snr', intensity_noise, n_channel_na)
%   uses the option 'na_snr' and the tuning parameter 'intensity_noise' to adjust the amplitude of the added Gaussian White noise. 
%   Paramter 'n_channel_na' determines the number of noise channels, and in default it is set equally with
%   the number of channels of the original signal,x. In related researches[1][2], 'intensity_noise'= 0 is recommended.    
% 
% imfs = namemd(x, ndir, stp_crit, stp_vec, 'na_fix', intensity_noise, n_channel_na)
%   uses the option 'na_fix' and the tuning parameter 'intensity_noise' to adjust the amplitude of the added Gaussian White noise. 
%   Paramter 'n_channel_na' determines the number of noise channels, and in default it is set equally with 
%   the number of channels of the original signal, x. In related researches[1][2], 'intensity_noise'= 1 is recommended.
% -----------------------
% This coding was also sucessfully applied in 4-channel-based EMG signal processing 
% (Zhang, Su, Xu, and Yao, Proc. 39th Annual Int. Conf. IEEE Eng.Medi.Bio.Soc.(IEEE EMBC), 2017; 
%  Zhang,Xu,Li,Duan,Wen,Yang,Zhang,and Yao,Biomedical Engineering Online, 2017)
%
% Copyright: Yuexin Wen, Yi Zhang, Peng Xu, Steven Su, Limei Xu, and Dezhong Yao (October 2017)
%
% Related researches:
% [1] Y. Zhang, P. Xu,P.Y. Li,K.Y.Duan,Y.X.Wen,Q.Yang,T.Zhang, and D.Z.Yao, Noise-assisted multivariate empirical mode decomposition
%     for multichannel EMG signals[J]. Biomedical Engineering Online, 2017, 16(1):107.
% [2] Y. Zhang, S. Su, P. Xu, D.Z Yao, Performance Evaluation of Noise-Assisted Multivariate Empirical Mode Decomposition and its Application to 
%    Multichannel EMG Signals, 39th Annual International Conference of the IEEE Engineering in Medicine and Biology Society, 2017, Jeju, Korea
% 
% Reference:
% [1]  Rehman and D. P. Mandic, "Multivariate Empirical Mode Decomposition", Proceedings of the Royal Society A, 2010
% [2]  G. Rilling, P. Flandrin and P. Goncalves, "On Empirical Mode Decomposition and its Algorithms", Proc of the IEEE-EURASIP
%      Workshop on Nonlinear Signal and Image Processing, NSIP-03, Grado (I), June 2003
% [3]  N. E. Huang et al., "A confidence limit for the Empirical Mode Decomposition and Hilbert spectral analysis",
%      Proceedings of the Royal Society A, Vol. 459, pp. 2317-2345, 2003
% -----------------------
% Demo:
%   x = randn(1000,3);
%   ndir = 64;
%   stp_crit = 'stop';
%   stp_vec = [0.075 0.75 0.075];
%   mode = 'na_fix';
%   intensity_noise = 1;
%   n_channel_na = 3;
%   imfs = namemd(x, ndir, stp_crit, stp_vec, mode, intensity_noise, n_channel_na);