function bandpass(OBJ)
%
%
%

% some pipelines get super fancy, but why get super fancy when you have an iir filter, because
% moar infinite

for i=1:length(OBJ)
	if OBJ(i).options.photometry.mod_bandpass
		for j=1:length(OBJ(i).traces)
			filter_data=OBJ(i).traces(j).raw;
			filter_data(isnan(filter_data))=0;
			[b,a]=ellip(5,.2,40,...
				[OBJ(i).metadata.traces(j).mod_freq-OBJ(i).options.photometry.mod_bandpass_bw/2 ...
				OBJ(i).metadata.traces(j).mod_freq+OBJ(i).options.photometry.mod_bandpass_bw/2]/(OBJ(i).metadata.fs/2),...
				'bandpass');
			OBJ(i).traces(j).raw_filt=filtfilt(b,a,filter_data);
		end
	end
end
