function [p]=localperiod(data,dt)

% The function PAZ generates a period using zero-crossing method 
% applied to data(n,k), where n specifies the length of time series, 
% and k is the number of IMFs.
% Non MATLAB Library routine used in the function is: FINDCRITICALPOINTS.
%
% Calling sequence-
% p=faz(data)
%
% Input-
%	data	- 2-D matrix of IMF components 
%	dt	    - time increment per point
% Output-
%	p	    - 2-D matrix f(n,k) that specifies frequency
%	
%
%----- Get dimensions
if nargin==1
    dt=1;
end

[f,a]=fazoi(data,dt);
p=1./f;