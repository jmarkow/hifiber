function demodulate(OBJ)
%
%
%

% use supplied references to demodulate each channel using lock-in principles

% (i.e. multiply by ref and 90 deg phase-shift)

% in other words, hypot it out and smooth yeahhhhhh

for i=1:length(OBJ)
	for j=1:length(OBJ(i).traces)

		% remove traces, and re-populate with the demodded version
		% do all to all, and bandpass as needed...

		for k=1:length(OBJ(i).references)


		end
	end
end
