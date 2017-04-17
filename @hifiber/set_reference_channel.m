function set_reference_channel(OBJ,CH,REF)
%
%
%

nchannels=length(OBJ.traces);

for i=1:length(OBJ)
	if (CH<=nchannels & CH>0) & (REF<=nchannels & REF>=0)
		OBJ.traces(CH).reference_channel=REF;
	end
end
