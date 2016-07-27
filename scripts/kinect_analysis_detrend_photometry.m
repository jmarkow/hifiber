

%% photometry analysis

gcamp_trace=mat(1,:);
gcamp_ts=mat(3,:);

win_size=5;
new_fs=100;
tau=5;
nidaq_fs=round(1./mean(diff(gcamp_ts)));

[b,a]=ellip(4,.2,40,[20]/(1e3/2),'low');
gcamp_trace_ds=downsample(filtfilt(b,a,-gcamp_trace),nidaq_fs/new_fs);
gcamp_ts_ds=downsample(gcamp_ts,10);
gcamp_trace_smooth=markolab_smooth(gcamp_trace_ds(:),round(tau*new_fs),'r','b');
%[b2,a2]=butter(3,[.01]/(new_fs/2),'low');
%gcamp_trace_smooth=filtfilt(b,a,gcamp_trace_ds(:));
%gcamp_trace_smooth=gcamp_trace_ds(:);
gcamp_trace_detrended=fluolab_detrend(gcamp_trace_smooth,'fs',new_fs,'win',win_size);
gcamp_trace_detrended=gcamp_trace_detrended(new_fs*win_size+tau*new_fs:end);
gcamp_ts_ds=gcamp_ts_ds(new_fs*win_size+tau*new_fs:end);


rcamp_trace=mat(2,:);

[b,a]=ellip(4,.2,40,[20]/(1e3/2),'low');
rcamp_trace_ds=downsample(filtfilt(b,a,-rcamp_trace),nidaq_fs/new_fs);
rcamp_trace_smooth=markolab_smooth(rcamp_trace_ds(:),round(tau*new_fs),'r','b');
%[b2,a2]=butter(3,[.01]/(new_fs/2),'low');
%gcamp_trace_smooth=filtfilt(b,a,gcamp_trace_ds(:));
%gcamp_trace_smooth=gcamp_trace_ds(:);
rcamp_trace_detrended=fluolab_detrend(rcamp_trace_smooth,'fs',new_fs,'win',win_size);
rcamp_trace_detrended=rcamp_trace_detrended(new_fs*win_size+tau*new_fs:end);
