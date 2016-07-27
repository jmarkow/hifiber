function [BOOTVALS_MAX,BOOTVALS_MIN,OBS_R,LAGS]=kinect_analysis_proc_photometry(DATA1,DATA2,varargin)
%
%
%
%
%

normalization='coeff';
nboots=1e3;
maxlags=100;

nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
	case 'normalization'
		normalization=varargin{i+1};
	case 'nboots'
		nboots=varargin{i+1};
	case 'maxlags'
		maxlags=varargin{i+1};
	otherwise
	end
end

BOOTVALS_MAX=nan(1,nboots);
BOOTVALS_MIN=nan(1,nboots);

for i=1:nboots
	scr=markolab_phase_scramble_1d(DATA1,0);
	r=xcorr(scr,DATA2,maxlags,normalization);
	BOOTVALS_MAX(i)=max(r);
	BOOTVALS_MIN(i)=min(r);
end

[OBS_R,LAGS]=xcorr(DATA1,DATA2,maxlags,normalization);
