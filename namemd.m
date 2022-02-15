function imfs = namemd(x, varargin)
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
% Copyright: Yuexin Wen, Yi Zhang, Steven Su, Peng Xu, and Dezhong Yao (October 2017)
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

% Initialization
[x, seq, ndir, stop_crit, stp_vec, stp_cnt, mode_memd, snr_add, std_add, N_channel_add] = ...
    set_value(x, nargin, varargin{:});

global N_length N_dim; 
N_dim_original = size(x, 2);

r = x;
t = 1:N_length;
n_imf = 1;
MAXITERATIONS = 200; 
MAX_N_IMFS = 100;
imfs_tmp = zeros(N_dim,MAX_N_IMFS,N_length);
imfs = cell(N_dim,1);

% Add noise
if strcmp(mode_memd,'na_snr');
    r = add_noise(x,mode_memd,snr_add,N_channel_add); 
    x = r;
elseif strcmp(mode_memd,'na_fix');
    r = add_noise(x,mode_memd,std_add,N_channel_add); 
    x = r;
end

% Decomposition 
while ~stop_emd(r, seq, ndir);
    m = r; % current mode   
    
    % computation of mean and stopping criterion
    if strcmp(stop_crit,'stop');
        [stop_sift,env_mean] = stop_sifting(m,t,stp_vec,seq,ndir);
    elseif strcmp(stop_crit,'fix_h'); 
        counter = 0;
        [stop_sift,env_mean,counter] = stop_sifting_fix(m,t,seq,ndir,stp_cnt,counter);
    end
    
    % sifting loop
    nbit = 1;
    while ~stop_sift && nbit<MAXITERATIONS;
        % sifting
        m = m - env_mean;
        % computation of mean and stopping criterion
        if strcmp(stop_crit,'stop');
            [stop_sift,env_mean] = stop_sifting(m,t,stp_vec,seq,ndir);
        elseif strcmp(stop_crit,'fix_h');
            [stop_sift,env_mean,counter] = stop_sifting_fix(m,t,seq,ndir,stp_cnt,counter);
        end
        nbit = nbit + 1;     
    end
    
    if nbit == MAXITERATIONS;
        warning('emd:warning','forced stop of sifting : too many iterations');
    end    
    
    imfs_tmp(:,n_imf,:) = m';
    n_imf = n_imf + 1;
    r = r - m;
    
    % In case the current mode is so small that machine precision can cause spurious extrema to appear
    if max(abs(r)) < (1e-3)*(max(abs(x)));
        if ~stop_sift;
            warning('emd:warning','forced stop of EMD : too small amplitude');
        else
            disp('forced stop of EMD : too small amplitude');
        end
        break
    end
end

% Return result 
imfs_tmp(:,n_imf,:) = r';
imfs_tmp = imfs_tmp(:,1:n_imf,:);
for i = 1 : N_dim;
    imfs{i} = reshape(imfs_tmp(i,:,:), size(imfs_tmp,2), size(imfs_tmp,3));
end
imfs = imfs(1:N_dim_original); 
end

%% Subfunctions

function stp = stop_emd(r, seq, ndir)

global N_dim;
ner = zeros(ndir,1);
for it=1:ndir
    if (N_dim~=3) % Multivariate signal (for N_dim ~=3) with hammersley sequence
        
        % Linear normalisation of hammersley sequence in the range of -1.00 ~ 1.00
        b = 2*seq(1:end,it) - 1;
        
        % Find angles corresponding to the normalised sequence
        tht = atan2(sqrt(flipud(cumsum(b(N_dim:-1:2).^2))),b(1:N_dim-1)).';
        % Find coordinates of unit direction vectors on n-sphere
        dir_vec(1:N_dim) = [1 cumprod(sin(tht))];
        dir_vec(1:N_dim-1) =  cos(tht) .*dir_vec(1:N_dim-1);

    else % 
         % Trivariate signal with hammersley sequence
         % Linear normalisation of hammersley sequence in the range of -1.0 ~ 1.0
        tt = 2*seq(1,it)-1;
        tt((tt>1)) = 1;
        tt((tt<-1)) = -1;
        
        % Normalize angle from 0 - 2*pi
        phirad = seq(2,it)*2*pi;
        st = sqrt(1.0-tt*tt);
        
        dir_vec(1) = st * cos(phirad);
        dir_vec(2) = st * sin(phirad);
        dir_vec(3) = tt;
    end
    % Projection of input signal on nth (out of total ndir) direction
    % vectors
    y = r * dir_vec';
    % Calculates the extrema of the projected signal
    [indmin, indmax] = local_peaks(y); % 
    
    ner(it) = length(indmin) + length(indmax); % 
end

% Stops if the all projected signals have less than 3 extrema
stp = all(ner < 3);% For vectors, ALL(V) returns logical 1 (TRUE) if none of the elements of the vector are zero. 
                   % For matrices, ALL(X) operates on the columns of X, returning a row vector of logical 1's and 0's. 
end

function [env_mean,nem,nzm,amp] = envelope_mean(m,t,seq,ndir)

global N_length N_dim;
NBSYM = 2;
count = 0;

env_mean = zeros(length(t),N_dim);
amp = zeros(length(t),1);
nem = zeros(ndir,1);
nzm = zeros(ndir,1);

for it = 1 : ndir; 
    if (N_dim ~=3) % Multivariate signal (for N_dim ~=3) with hammersley sequence
        % Linear normalisation of hammersley sequence in the range of -1.00 - 1.00
        b = 2*seq(1:end,it)-1;
        % Find angles corresponding to the normalised sequence
        tht = atan2(sqrt(flipud(cumsum(b(N_dim:-1:2).^2))),b(1:N_dim-1)).';
        % Find coordinates of unit direction vectors on n-sphere
        dir_vec(1:N_dim) = [1 cumprod(sin(tht))];
        dir_vec(1:N_dim-1) =  cos(tht).*dir_vec(1:N_dim-1);
    else % Trivariate signal with hammersley sequence£¬
        % Linear normalisation of hammersley sequence in the range of -1.0 - 1.0
        tt = 2*seq(1,it)-1;
        tt((tt>1)) = 1;
        tt((tt<-1)) = -1;
        
        % Normalize angle from 0 - 2*pi
        phirad = seq(2,it)*2*pi;
        st = sqrt(1.0-tt*tt);
        
        dir_vec(1)=st * cos(phirad);
        dir_vec(2)=st * sin(phirad);
        dir_vec(3)=tt;
    end
    
    % Projection of input signal on nth (out of total ndir) direction vectors
    y(1:N_length) = dir_vec * m(1:N_length,:)';
    
    % Calculates the extrema of the projected signal
    [indmin, indmax] = local_peaks(y);
    
    nem(it) = length(indmin) + length(indmax); 
    
    indzer = zero_crossings(y); 
    nzm(it) = length(indzer); 
    
   
    [tmin,tmax,zmin,zmax,mode] = boundary_conditions(indmin,indmax,t,y,m,NBSYM);
    
    % Calculate multidimensional envelopes using spline interpolation
    % Only done if number of extrema of the projected signal exceed 3
    if(mode)
        env_min = spline(tmin,zmin.',t).'; % lower envelopes
        env_max = spline(tmax,zmax.',t).'; % upper envelopes 
        
        amp = amp + sqrt(sum((env_max-env_min).^2,2))/2; 
        env_mean = env_mean + (env_max+env_min)/2; 
    else % if the projected signal has inadequate extrema
        count = count + 1; 
    end
end

if ndir > count; 
    env_mean = env_mean/(ndir-count); 
    amp = amp/(ndir-count);
else 
    env_mean = zeros(N_length,N_dim);
    amp = zeros(N_length,1);
    nem = zeros(1,ndir);
end
end

function [stp,env_mean,nzm] = stop_sifting(m,t,stp_vec,seq,ndir)

global N_length N_dim;
sd = stp_vec(1);
sd2 = stp_vec(2);
tol = stp_vec(3);
   
try
    [env_mean, nem, nzm, amp] = envelope_mean(m,t,seq,ndir);
    sx = sqrt(sum(env_mean.^2,2));
    if(amp) % something is wrong here
        sx = sx./amp; % amp = amp/(ndir-count);
    end
    stp = ~( (mean(sx>sd) > tol || any(sx > sd2)) && any(nem > 2) );
catch
    env_mean = zeros(N_length,N_dim);
    stp = 1;
end
end

function [stp,env_mean,counter]= stop_sifting_fix(m,t,seq,ndir,stp_count,counter)

global N_length N_dim;
try
    [env_mean,nem,nzm] = envelope_mean(m,t,seq,ndir);
    if (all(abs(nzm-nem)>1)) 
        stp = 0;
        counter = 0;
    else
        counter = counter+1;
        stp = (counter >= stp_count);
    end
catch
    env_mean = zeros(N_length,N_dim);
    stp = 1;
end
end

function [tmin,tmax,zmin,zmax,mode] = boundary_conditions(indmin,indmax,t,x,z,nbsym)

lx = length(x);
if (length(indmin) + length(indmax) < 3)
    mode = 0;
    tmin=NaN;tmax=NaN;zmin=NaN;zmax=NaN;
    return
else
    mode = 1; %the projected signal has inadequate extrema
end

% boundary conditions for interpolations :
if indmax(1) < indmin(1)
    if x(1) > x(indmin(1))
        lmax = fliplr(indmax(2:min(end,nbsym+1))); % Flip matrix left to right
        lmin = fliplr(indmin(1:min(end,nbsym))); 
        lsym = indmax(1);
    else
        lmax = fliplr(indmax(1:min(end,nbsym)));
        lmin = [fliplr(indmin(1:min(end,nbsym-1))),1];
        lsym = 1;
    end
else
    if x(1) < x(indmax(1))
        lmax = fliplr(indmax(1:min(end,nbsym)));
        lmin = fliplr(indmin(2:min(end,nbsym+1)));
        lsym = indmin(1);
    else
        lmax = [fliplr(indmax(1:min(end,nbsym-1))),1];
        lmin = fliplr(indmin(1:min(end,nbsym)));
        lsym = 1;
    end
end

if indmax(end) < indmin(end)
    if x(end) < x(indmax(end))
        rmax = fliplr(indmax(max(end-nbsym+1,1):end));
        rmin = fliplr(indmin(max(end-nbsym,1):end-1));
        rsym = indmin(end);
    else
        rmax = [lx,fliplr(indmax(max(end-nbsym+2,1):end))];
        rmin = fliplr(indmin(max(end-nbsym+1,1):end));
        rsym = lx;
    end
else
    if x(end) > x(indmin(end))
        rmax = fliplr(indmax(max(end-nbsym,1):end-1));
        rmin = fliplr(indmin(max(end-nbsym+1,1):end));
        rsym = indmax(end);
    else
        rmax = fliplr(indmax(max(end-nbsym+1,1):end));
        rmin = [lx,fliplr(indmin(max(end-nbsym+2,1):end))];
        rsym = lx;
    end
end
tlmin = 2*t(lsym)-t(lmin);
tlmax = 2*t(lsym)-t(lmax);
trmin = 2*t(rsym)-t(rmin);
trmax = 2*t(rsym)-t(rmax);

% in case symmetrized parts do not extend enough
if tlmin(1) > t(1) || tlmax(1) > t(1)
    if lsym == indmax(1)
        lmax = fliplr(indmax(1:min(end,nbsym)));
    else
        lmin = fliplr(indmin(1:min(end,nbsym)));
    end
    if lsym == 1
        error('bug')
    end
    lsym = 1;
    tlmin = 2*t(lsym)-t(lmin);
    tlmax = 2*t(lsym)-t(lmax);
end

if trmin(end) < t(lx) || trmax(end) < t(lx)
    if rsym == indmax(end)
        rmax = fliplr(indmax(max(end-nbsym+1,1):end));
    else
        rmin = fliplr(indmin(max(end-nbsym+1,1):end));
    end
    if rsym == lx
        error('bug')
    end
    rsym = lx;
    trmin = 2*t(rsym)-t(rmin);
    trmax = 2*t(rsym)-t(rmax);
end
zlmax =z(lmax,:);
zlmin =z(lmin,:);
zrmax =z(rmax,:);
zrmin =z(rmin,:);

tmin = [tlmin t(indmin) trmin];
tmax = [tlmax t(indmax) trmax];
zmin = [zlmin; z(indmin,:); zrmin];
zmax = [zlmax; z(indmax,:); zrmax];
end

function [indmin, indmax] = local_peaks(x)

if(all(x < 1e-5))
    x=zeros(1,length(x));
end
m = length(x);

% Calculates the extrema of the projected signal
% Difference between subsequent elements:
dy = diff(x); % difference of signal 
a = find(dy~=0); 
lm = find(diff(a)~=1) + 1; 
d = a(lm) - a(lm-1); 
a(lm) = a(lm) - floor(d/2);
a(end+1) = m;
ya  = x(a); 

if(length(ya) > 1)
    % Maxima
    [pks_max,loc_max] = peaks(ya);
    % Minima
    [pks_min,loc_min] = peaks(-ya);
    
    if(~isempty(pks_min))
        indmin = a(loc_min);
    else
        indmin = NaN;
    end
    
    if(~isempty(pks_max))
        indmax = a(loc_max);
    else
        indmax = NaN;
    end
else
    indmin=NaN;
    indmax=NaN;
end
end

function [pks_max,locs_max] = peaks(X)
% Find the maximum and where it is located

dX = sign(diff(X));
locs_max = find((dX(1:end-1) >0) & (dX(2:end) <0)) + 1;
pks_max = X(locs_max);
end

function indzer = zero_crossings(x)
indzer = find(x(1:end-1).*x(2:end)<0);
if any(x == 0)
    iz = find( x==0 );
    if any(diff(iz)==1)
        zer = x == 0;
        dz = diff([0 zer 0]);
        debz = find(dz == 1);
        finz = find(dz == -1)-1;
        indz = round((debz+finz)/2);
    else
        indz = iz;
    end
    indzer = sort([indzer indz]);
end
end

function seq = hamm(n,base)
% Generate Hammersley sequence 

seq = zeros(1,n);
if ( 1 < base )
    seed = 1:1:n;
    base_inv = inv(base);
    while ( any ( seed ~= 0 ) )
        digit = mod (seed(1:n), base);
        seq = seq + digit * base_inv;
        base_inv = base_inv / base;
        seed = floor (seed / base );
    end
else
    temp = 1:1:n;
    seq = (mod(temp,(-base + 1))+0.5)/(-base);
end
end

function [x,seq,ndir,stp_crit,stp_vec,stp_cnt,mode_memd,snr_add,std_add,N_channel_add] = set_value(x,N_nargin,varargin)
global N_length N_dim; 
ndir = [];
stp_crit = [];
stp_vec = [];
stp_cnt  = [];
mode_memd = [];
std_add = []; 
snr_add = [];

if N_nargin < 1 || N_nargin > 7;
    InputError_Number = generatemsgid('InputError_Number');
    error(InputError_Number, 'number of input arguments is invalid');
end

%%%%%%%%%%
x = double(x);

% Rescale input signal if required
if (any(size(x)) == 0)
    datamsgid = generatemsgid('emptyDataSet');
    error(datamsgid,'Data set cannot be empty.');
end

%%%%%%%%%%
N_dim = size(x,2); 

% Specifies maximum number of channels that can be processed by the code
% Its maximum possible value is 32
Max_channels = 5100; %default value was 16 (Changed by Hamed)

if(N_dim < 2 || N_dim > Max_channels)
    error('Function only processes the signal having 2 and 16 channels.');
end

% Length of data
N_length = size(x,1);

%%%%%%%%%%
if N_nargin >= 2;
    ndir = varargin{1}; % Projective dimension of signal 
end

if N_nargin >= 3;
    stp_crit = varargin{2};
end

if N_nargin >= 4;
    if strcmp(varargin{2},'fix_h'); 
        stp_cnt = varargin{3};
    elseif strcmp(varargin{2},'stop');
        stp_vec = varargin{3};
    elseif ~xor(strcmp(varargin{2},'fix_h'),strcmp(varargin{2},'stop'));
        Nmsgid = generatemsgid('invalid stop_criteria');
        error(Nmsgid,'stop_criteria should be either fix_h or stop');
    end
end

if N_nargin >= 5 && ischar(varargin{4});
    mode_memd = varargin{4};
end

if ~ischar(mode_memd);
    mode_memd = 'memd'; % default
elseif ~strcmp(mode_memd,'memd') && ~strcmp(mode_memd,'na_snr') && ~strcmp(mode_memd,'na_fix');
    mode_memd = 'memd'; % default
end

%%%%%%%%%%% Noise arguments 
if N_nargin >= 6 && isnumeric(varargin{5});
    if strcmp(mode_memd,'na_snr'); 
        snr_add = varargin{5};
    elseif strcmp(mode_memd,'na_fix');
        std_add = varargin{5};
    end
end

if strcmp(mode_memd,'na_fix') || strcmp(mode_memd,'na_snr');
    if N_nargin >= 7 && isnumeric(varargin{6}) && ~isempty(varargin{6});
        N_channel_add = varargin{6};
    else
        N_channel_add = N_dim;
    end
    N_dim = N_dim + N_channel_add; % Update the dimension of data
else
    N_channel_add = [];
end

%%%%%%%%%% Validation of the input arguments 
if ~isempty(ndir) && (~isnumeric(ndir) || ~isscalar(ndir) || any(rem(ndir,1)) || (ndir < 2))
    Nmsgid = generatemsgid('invalid num_dir');
    error(Nmsgid,'num_dir should be an integer greater than or equal to 6.');
end

if ~isempty(stp_crit) && (~ischar(stp_crit) || ~xor(strcmp(stp_crit,'fix_h'),strcmp(stp_crit,'stop')))
    Nmsgid = generatemsgid('invalid stop_criteria');
    error(Nmsgid,'stop_criteria should be either fix_h or stop');
end

if ~isempty(stp_vec) && (~isnumeric(stp_vec) || length(stp_vec)~=3 || ~strcmp(stp_crit,'stop'))
    Nmsgid = generatemsgid('invalid stop_vector');
    error(Nmsgid,'stop_vector should be an array with three elements e.g. default is [0.075 0.75 0.075] ');
end

if ~isempty(stp_cnt) && (~isnumeric(stp_cnt) || ~isscalar(stp_cnt) || ...
        any(rem(stp_cnt,1)) || (stp_cnt < 0) || ~strcmp(stp_crit,'fix_h'));
    Nmsgid = generatemsgid('invalid stop_count');
    error(Nmsgid,'stop_count should be a nonnegative integer');
end

%%%%%%%%%% Default values
if isempty(ndir);
    ndir = 64; % default
end

if isempty(stp_crit);
    stp_crit = 'stop'; % default
end

if isempty(stp_vec);
    stp_vec = [0.075,0.75,0.075]; % default
end

if isempty(stp_cnt);
    stp_cnt = 2; % default
end

%%%%%%%%%% Initializations for Hammersley function
%prm = [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,...
%    101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199, ...
%    211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,293,307,311,313,317,331]; 
prm = primes(90000); % modified by hamed
% Find the pointset for the given input signal
if (N_dim==3)
    base(1) = -ndir;
    base(2) = 2;
    seq = zeros(N_dim-1,ndir);
    for it = 1 : N_dim-1;
        seq(it,:) = hamm(ndir,base(it));
    end
else
    base = zeros(1,N_dim);
    base(1) = -ndir;
    for iter = 2 : N_dim;
        base(iter) = prm(iter-1);
    end
    seq = zeros(N_dim,ndir);
    for it = 1 : N_dim;
        seq(it,:) = hamm(ndir,base(it));
    end
end
end
