function rereference(OBJ)
% reference all channels that haven't yet been rereferenced
%

% check the metadata for channels that have been assigned references

to_reref=find(cellfun(@(x) ~isempty(x),{OBJ.metadata.channels(:).reference_channel}));

for i=to_reref

	switch lower(OBJ.options.photometry.rereference_method(1))

	case 'v'

		% vector rejection
		% make sure that the baseline has been subtracted first (slow drift should
		% not factor into this calculation)


	case 'l'

		% least squares (should give same answer, leaving here for completeness)

	otherwise
		error('Did not understand rereferencing method ((v)ector and (l)east squares)');
	end

end
