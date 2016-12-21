function set_reference_channel(OBJ,CH,REF)
%
%
%

nchannels=numel(OBJ.metadata.channels);

if (CH<=nchannels & CH>0) & (REF<=nchannels & REF>=0)
	OBJ.metadata.channels(CH).reference_channel=REF;
end
