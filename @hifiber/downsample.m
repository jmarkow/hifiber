function downsample(OBJ)
% downsamples data for downstream processing
%

% anti-alias filtering

valid_types={'dff','raw','baseline','baseline_rem','reref'};

fprintf('Downsampling data...\n');

for i=1:length(OBJ)

	% if mod(OBJ(i).metadata.fs,OBJ(i).options.new_fs)~=0
	% 	error('Downsampling is only supported for integer downsampling factors');
	% end

	% [b,a]=ellip(4,.2,40,[.75*(OBJ(i).options.new_fs/2)]/(OBJ(i).metadata.fs/2),'low');
	% downsample_factor=OBJ(i).metadata.fs/OBJ(i).options.new_fs;
	upd=hifiber.proc_timer(length(OBJ(i).traces));

	[p,q]=rat(OBJ(i).options.new_fs,OBJ(i).metadata.fs);
	nans=isnan(OBJ(i).timestamps);

	for j=1:length(OBJ(i).traces)

		data_types=fieldnames(OBJ(i).traces(j));

		for k=1:length(data_types)

			if any(strcmp(data_types{k},valid_types))

			% let MATLAB/god sort 'em out

				%OBJ(i).traces(j).(data_types{k})(nans)=nan;
				signal_len=length(OBJ(i).traces(j).(data_types{k}));
				[OBJ(i).traces(j).(data_types{k}),new_ts,b]=resample(OBJ(i).traces(j).(data_types{k}),...
					[0:signal_len-1]/OBJ(i).metadata.fs,OBJ(i).options.new_fs);

				% guess what punk, some more *edge effects since the signal is assumed to go to zero at the edges, nan out
				% taps/2 or pad intelligently

				filt_length=numel(b)/OBJ(i).metadata.fs;
				filt_length_res=round(filt_length*OBJ(i).options.new_fs)+10;

				OBJ(i).traces(j).(data_types{k})(1:filt_length_res)=nan;
				OBJ(i).traces(j).(data_types{k})(end-filt_length_res:end)=nan;

			end
		end

		upd(j);

	end

	new_timestamps=OBJ(i).timestamps(round(new_ts*OBJ(i).metadata.fs)+1);
	OBJ(i).timestamps=new_timestamps;

	%OBJ(i).timestamps=downsample(OBJ(i).timestamps,downsample_factor);
	OBJ(i).metadata.fs=OBJ.options.new_fs;

end
