function detect_references(OBJ)
% each channel should have mod freq associated with it, we then
% take care of the rest...be our guest
%
%
%

for i=1:length(OBJ)
	for j=1:length(OBJ(i).traces)

		% think about what kind of data structures we need for this crap
		% maybe stash the reference signal so we can double check things post-hoc

		% bandpass about the mod freq, nice 'n sleazy

		% put that in here rather than bandpassing first, otherwise cross-referencing gets confusing

		use_data=OBJ(i).traces(j).raw;
		nans=isnan(use_data);

		if OBJ(i).options.photometry.mod_bandpass
			use_data(nans)=0;
			[b,a]=ellip(5,.2,40,...
				[OBJ(i).metadata.traces(j).mod_freq-OBJ(i).options.photometry.mod_bandpass_bw/2 ...
				OBJ(i).metadata.traces(j).mod_freq+OBJ(i).options.photometry.mod_bandpass_bw/2]/(OBJ(i).metadata.fs/2),...
				'bandpass');
			use_data=filtfilt(b,a,use_data);
			use_data(nans)=nan;
		end

		tvec=[0:numel(OBJ(i).traces(j).raw)-1]/OBJ(i).metadata.fs;

		[params,fit_fun]=photometry.get_demod_reference(use_data,...
			tvec,OBJ(i).metadata.fs,OBJ(i).metadata.traces(j).mod_freq);

		% reset amplitude to 1

		OBJ(i).references(j).x=fit_fun([1 params(2:end)],tvec);
		OBJ(i).references(j).y=fit_fun([1 params(2) params(3)+pi/2 params(4)],tvec);

	end
end
