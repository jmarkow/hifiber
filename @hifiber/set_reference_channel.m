function set_reference_channel(OBJ,CH,REF)
%
%
%

nchannels=length(OBJ.traces);

if (CH<=nchannels & CH>0) & (REF<=nchannels & REF>=0)
	OBJ.traces(CH).reference_channel=REF;
end
