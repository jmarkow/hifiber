function subtract_baseline(OBJ)
%
%

for i=1:length(OBJ.traces)
	OBJ.traces(i).baseline_rem=OBJ.traces(i).raw-OBJ.traces(i).baseline;
end
