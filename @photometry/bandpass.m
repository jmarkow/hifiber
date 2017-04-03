function bandpass(OBJ)
%
%
%

% some pipelines get super fancy, but why get super fancy when you have an iir filter, because
% moar infinite

for i=1:length(OBJ)
	if OBJ(i).options.mod_bandpass
		[b,a]=ellip(5,.2,40,...
			[OBJ(i).traces(j).mod_freq-OBJ(i).options.mod_bandpass/2 ...
			OBJ(i).traces(j).mod_freq+OBJ(i).options.mod_bandpass/2]/(OBJ(i).metadata.fs/2),...
			'bandpass');
		OBJ(i).traces(j).raw_filt=filtfilt(b,a,OBJ(i).traces(j).raw);
	end
end
