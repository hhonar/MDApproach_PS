% Subroutine for creating grided covariance subplots
function [h] = gsplot(A,varargin)
[X,Y] = meshgrid(0:1:size(A,1), 0:1:size(A,2));
h = surface(X,-Y,padarray(A,[1 1],'post')','LineStyle','-');

p = inputParser;
addRequired(p,'A',@ismatrix);

defaultcolor = [1 1 1];   % add the optional grid color with this default value
addOptional(p,'EdgeColor',defaultcolor,@isvector)

defaultLW = 0.0001;   % add the optional window size  with this default value
addParameter(p,'LineWidth',defaultLW,@isnumeric)

parse(p,A,varargin{:});
A   = p.Results.A;
color = p.Results.EdgeColor;
LW = p.Results.LineWidth;
%option  = p.Results.option;
%wintype = p.Results.wintype;

set(h,'EdgeColor',color);
set(h,'LineWidth',LW);
axis square;%colorbar;caxis([-1 1]);
set(gca,'visible','off')
%set(h, 'XTick', [], 'YTick', [], 'CLim', [-1 1], 'XColor', [1 1 1], 'YColor', [1 1 1]);
return;
end