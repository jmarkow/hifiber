function downsample(OBJ)
% downsamples data for downstream processing
%

% anti-alias filtering


fprintf('Downsampling data...\n');

for i=1:length(OBJ)

	if mod(OBJ(i).metadata.fs,OBJ(i).options.photometry.new_fs)~=0
		error('Downsampling is only supported for integer downsampling factors');
	end

	[b,a]=ellip(4,.2,40,[.75*(OBJ(i).options.photometry.new_fs/2)]/(OBJ(i).metadata.fs/2),'low');
	downsample_factor=OBJ(i).metadata.fs/OBJ(i).options.photometry.new_fs;
	upd=photometry.proc_timer(length(OBJ(i).traces));

	for j=1:length(OBJ(i).traces)
		data_types=fieldnames(OBJ(i).traces(j));
		for k=1:length(data_types)
			OBJ(i).traces(j).(data_types{k})=filtfilt(b,a,OBJ(i).traces(j).(data_types{k}));
			OBJ(i).traces(j).(data_types{k})=downsample(OBJ(i).traces(j).(data_types{k}),downsample_factor);
		end
		upd(j);
	end

	OBJ(i).timestamps=downsample(OBJ(i).timestamps,downsample_factor);
	OBJ(i).metadata.fs=OBJ.options.photometry.new_fs;

end
