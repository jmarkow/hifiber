function set_channel_names(OBJ,NAMES)
%
%
%

if iscell(NAMES) & length(NAMES)==length(OBJ.metadata.channels)
	for i=1:length(NAMES)
		OBJ.metadata.channels(i).name=NAMES{i};
	end
end
