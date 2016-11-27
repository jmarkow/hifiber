%% photometry analysis from Minsuk's GCaMP/RCaMP data

win_size=2;
new_fs=100;
tau=1;
nidaq_fs=1e3;

traces=[AD2_12.data(:) AD3_12.data(:)];
ts=[1:size(traces,1)]/nidaq_fs;

[b,a]=ellip(4,.2,40,[20]/(1e3/2),'low');
traces_ds=downsample(filtfilt(b,a,traces),nidaq_fs/new_fs);
ts_ds=downsample(ts,nidaq_fs/new_fs);
traces_smooth=markolab_smooth(traces_ds,round(tau*new_fs),'r','e');
%[b2,a2]=butter(3,[.1]/(new_fs/2),'low');
%traces_smooth=filtfilt(b2,a2,traces_ds);
%gcamp_trace_smooth=gcamp_trace_ds(:);

traces_detrended=[];
traces_detrended(:,1)=fluolab_detrend(traces_smooth(:,1),'fs',new_fs,'win',win_size);
traces_detrended(:,2)=fluolab_detrend(traces_smooth(:,2),'fs',new_fs,'win',win_size);

traces_detrended=traces_detrended(new_fs*win_size:end,:);
ts_ds=ts_ds(new_fs*win_size:end);

%%

figure();plot(ts_ds,traces_detrended(:,1),'g-');
hold on
plot(ts_ds,traces_detrended(:,2),'r-');
