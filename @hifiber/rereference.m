function rereference(OBJ)
% reference all channels that haven't yet been rereferenced
%

% check the metadata for channels that have been assigned references

to_reref=find(cellfun(@(x) ~isempty(x),{OBJ.traces(:).reference_channel}));

for i=to_reref

	switch lower(OBJ.options.rereference_method(1))

	case 'v'

		% vector rejection
		% make sure that the baseline has been subtracted first (slow drift should
		% not factor into this calculation)

		sig=OBJ.traces(i).baseline_rem;
		reference=OBJ.traces(OBJ.traces(i).reference_channel).baseline_rem;

		use_samples=~(isnan(sig)|isnan(reference));

		num=sig(use_samples)'*reference(use_samples);
		den=reference(use_samples)'*reference(use_samples);
		OBJ.traces(i).reref=sig-reference*(num/den);

	case 'l'

		% least squares (should give same answer, leaving here for completeness)

	otherwise
		error('Did not understand rereferencing method ((v)ector and (l)east squares)');
	end

end
