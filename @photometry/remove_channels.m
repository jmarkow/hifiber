function remove_channels(OBJ,NAMES)
%
%
%
%

to_rem=false(1,length(OBJ.metadata.channels));

for i=1:length(NAMES)
	idx=cellfun(@(x) ~isempty(x),strfind({OBJ.metadata.channels(:).name},NAMES{i}));
	to_rem(idx)=true;
end

OBJ.traces(to_rem)=[];
OBJ.metadata.channels(to_rem)=[];

% to_rem
