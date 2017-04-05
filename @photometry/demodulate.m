function demodulate(OBJ)
%
%
%

% use supplied references to demodulate each channel using lock-in principles

% (i.e. multiply by ref and 90 deg phase-shift)

% in other words, hypot it out and smooth yeahhhhhh

for i=1:length(OBJ)

	demod_samples=round(OBJ(i).metadata.fs*OBJ(i).options.photometry.demod_tau);
	demod_kernel=ones(demod_samples,1)/demod_samples;

	traces=cat(2,OBJ(i).traces(:).raw);
	counter=1;

	for j=1:length(OBJ(i).traces)

		% remove traces, and re-populate with the demodded version
		% do all to all, and bandpass as needed...

		for k=1:length(OBJ(i).references)

			% bandpass for the reference fs

			use_data=photometry.bandpass(traces(:,j),OBJ(i).metadata.traces(k).mod_freq,...
				OBJ(i).options.photometry.mod_bandpass_bw,OBJ(i).metadata.fs);
			OBJ(i).traces(counter).raw=conv(hypot(use_data.*OBJ(i).references(k).x,use_data.*OBJ(i).references(k).y),demod_kernel,'same');
			counter=counter+1;

		end
	end
end
