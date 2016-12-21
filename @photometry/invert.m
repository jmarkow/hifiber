function invert(OBJ)
%
%
%

for i=1:length(OBJ.traces)
	OBJ.traces(i).raw=-OBJ.traces(i).raw;
end
