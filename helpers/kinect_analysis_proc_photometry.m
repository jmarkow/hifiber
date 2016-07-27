function [PROC_DATA,PROC_TS]=kinect_analysis_proc_photometry(DATA,TS,varargin)
%
%
%
%
%

smooth_type='e';
smooth_tau=.3;
new_fs=100;
detrend_win_size=3;
dff=1;
detrend=true;
nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
  switch lower(varargin{i})
    case 'smooth_type'
      smooth_type=varargin{i+1};
    case 'smooth_tau'
      smooth_tau=varargin{i+1};
    case 'new_fs'
      new_fs=varargin{i+1};
    case 'dff'
      dff=varargin{i+1};
    otherwise
  end
end

fs=round(1./mean(diff(TS)));
new_fs=100;
down_fact=round(fs/new_fs);

[b,a]=ellip(4,.2,40,[20]/(fs/2),'low');
gcamp_trace_ds=downsample(filtfilt(b,a,DATA),down_fact);
PROC_TS=downsample(TS,down_fact);

if smooth_tau>0
    gcamp_trace_smooth=markolab_smooth(gcamp_trace_ds(:),round(smooth_tau*new_fs),'r',smooth_type);
else
    gcamp_trace_smooth=gcamp_trace_ds(:);
end

if detrend
    PROC_DATA=fluolab_detrend(gcamp_trace_smooth,'fs',new_fs,'win',detrend_win_size,'dff',dff);
else
    PROC_DATA=gcamp_trace_smooth;
end
PROC_DATA=PROC_DATA(new_fs*detrend_win_size:end);
PROC_TS=PROC_TS(new_fs*detrend_win_size:end);
