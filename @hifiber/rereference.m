function rereference(OBJ)
% reference all channels that haven't yet been rereferenced
%

% check the metadata for channels that have been assigned references

for i=1:length(OBJ)
	to_reref=find(cellfun(@(x) ~isempty(x),{OBJ(i).traces(:).reference_channel}));
	for j=to_reref

		switch lower(OBJ(i).options.rereference_method(1))

			case 'v'

			% vector rejection
			% make sure that the baseline has been subtracted first (slow drift should
			% not factor into this calculation)

			sig=OBJ(i).traces(j).baseline_rem;
			reference=OBJ(i).traces(OBJ(i).traces(j).reference_channel).baseline_rem;

			use_samples=~(isnan(sig)|isnan(reference));

			num=sig(use_samples)'*reference(use_samples);
			den=reference(use_samples)'*reference(use_samples);
			OBJ(i).traces(j).reref=sig-reference*(num/den);
			OBJ(i).traces(j).dff=OBJ(i).traces(j).reref./OBJ(i).traces(j).baseline;


			case 'l'

			% least squares (should give same answer, leaving here for completeness)

			otherwise
			error('Did not understand rereferencing method ((v)ector and (l)east squares)');
		end

	end
end
