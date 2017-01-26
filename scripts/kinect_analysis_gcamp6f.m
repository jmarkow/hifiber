%%% script for analyzing the gcamp6f mice

% first use date cutoffs for when we did raw and modulated recordings, processing should be the same...
% 2016-09-02, everything after this two channels, first green second red, inverted during acquisition
% 2016-08-02, raw 4 channels, 1 and 3 green and red respectively
% 2016-07-19, raw 4 channels, 1 and 3 green and red respectively
% 2016-06-07, dc offset scale 100, 4 channels

% consider re-combining dc offset...simple enough to process everything after 2016-07-19 for now...


times=zeros(1,length(extract_object));

for i=1:length(extract_object)
	times(i)=datenum(extract_object(i).metadata.extract.StartTime,'yyyy-mm-dd');

end

cutoffs=[ datenum('2016-06-07') datenum('2016-07-19','yyyy-mm-dd') datenum('2016-09-02','yyyy-mm-dd') inf ];
[~,phot_category]=histc(times,cutoffs);

% bin 1 we need to reincoporate the dc offset
% bin 2 just take channels 1 and 3
% bin 3 same, but take channels 1 and 2

for i=find(phot_category==2)

	extract_object(i).neural_data.photometry.set_channel_names(extract_object(i).metadata.extract.NidaqChannelNames);
	extract_object(i).neural_data.photometry.remove_channels({'ai1','ai3'});

		% if length(extract_object(i).neural_data.photometry.traces)>2
		% 	extract_object(i).neural_data.photometry.traces([2 4])=[];
		% 	extract_object(i).neural_data.photometry.metadata.channels([2 4])=[];
		% end
		%
		% % check for baseline, baseline subtr and dff
		%

	extract_object(i).neural_data.photometry.get_baseline;
	extract_object(i).neural_data.photometry.subtract_baseline;
	extract_object(i).neural_data.photometry.get_dff;

end