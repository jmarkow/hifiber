function set_modulation_frequency(OBJ,FREQ,CH)
%
%
%
%

for i=1:length(OBJ)
	if CH<=length(OBJ(i).traces)
		OBJ(i).traces(CH).mod_freq=FREQ;
	end
end
