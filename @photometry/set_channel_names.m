function set_channel_names(OBJ,NAMES)
%
%
%

if iscell(NAMES) & length(NAMES)==length(OBJ.metadata.traces)
	for i=1:length(NAMES)
		OBJ.metadata.traces(i).name=NAMES{i};
	end
end
