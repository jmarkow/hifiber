function [RPS,SCORE,CHANGEPOINTS,THRESH]=kinect_analysis_changepoints(FEATURES,varargin)
%
%
%
%

if ndims(squeeze(FEATURES))==3
	[r,c,z]=size(FEATURES);
	FEATURES=reshape(FEATURES,r*c,z)';
end

delta_thresh=.3;
delta_win=4;
smooth_sig=2;
jl_bound_eps=.25;

nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
  switch lower(varargin{i})
    case 'delta_thresh'
      delta_thresh=varargin{i+1};
    case 'delta_win'
      delta_win=varargin{i+1};
    case 'smooth_sig'
      smooth_sig=varargin{i+1};
    otherwise
  end
end

n_components=round(4*log(size(FEATURES,2))/(jl_bound_eps^2/2-jl_bound_eps^3/3));
%n_components=400;
RPS=kinect_gaussproj(zscore(FEATURES),n_components);
RPS=zscore(RPS');

RPS_deltas=markolab_deltacoef(RPS,delta_win); % delta coefficients, lag of 4 frames
SCORE=sum(abs(RPS_deltas)>delta_thresh); % binarize deltas
h=normpdf([-10:10],0,smooth_sig); % gauss smooth
SCORE=conv(SCORE,h,'same'); % convolve
THRESH=mean(SCORE)+.5*std(SCORE);
idx=1:length(SCORE)-1;

[~,CHANGEPOINTS]=findpeaks(SCORE,'minpeakheight',THRESH);
