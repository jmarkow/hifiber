function set_channel_names(OBJ,NAMES,NUM)
%
%
%

if nargin<3
	NUM=[];
end

if iscell(NAMES) & length(NAMES)==length(OBJ.metadata.traces)
	for i=1:length(NAMES)
		OBJ.metadata.traces(i).name=NAMES{i};
	end
elseif isa(NAMES,'string') & ~isempty(NUM)
	OBJ.metadata.traces(NUM).name=NAMES;
end
