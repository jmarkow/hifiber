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

		tvec=[0:numel(OBJ(i).traces(j).raw-1]/OBJ(i).metadata.fs;

		if isfield(OBJ(i).traces(j),'raw_filt')
			use_field='raw_filt';
		else
			use_field='raw';
		end

		[params,fit_fun]=photometry.get_demod_reference(OBJ(i).traces(j).(use_field),tvec,OBJ(i).traces(j).mod_freq);

		% reset amplitude to 1

		OBJ(i).references.x(j)=fit_fun([1 params(2:end)],tvec);
		OBJ(i).references.y(j)=fit_fun([1 params(2) params(3)+pi/2 params(4)],tvec);

	end
end
