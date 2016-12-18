function downsample(OBJ)
% downsamples data for downstream processing
%

% anti-alias filtering

if mod(OBJ.metadata.fs,OBJ.options.photometry.new_fs)~=0
	error('Downsampling is only supported for integer downsampling factors');
end

fprintf('Downsampling data...\n');

[b,a]=ellip(4,.2,40,[.75*(OBJ.options.photometry.new_fs/2)]/(OBJ.metadata.fs/2),'low');
downsample_factor=OBJ.metadata.fs/OBJ.options.photometry.new_fs;
upd=photometry.proc_timer(length(OBJ.traces));

for i=1:length(OBJ.traces)
	data_types=fieldnames(OBJ.traces(i));
	for j=1:length(data_types)
		OBJ.traces(i).(data_types{j})=filtfilt(b,a,OBJ.traces(i).(data_types{j}));
		OBJ.traces(i).(data_types{j})=downsample(OBJ.traces(i).(data_types{j}),downsample_factor);
	end
	upd(i);
end

OBJ.timestamps=downsample(OBJ.timestamps,downsample_factor);
OBJ.metadata.fs=OBJ.options.photometry.new_fs;
