function FILT_DATA=bandpass(DATA,FC,BW,FS)
%
%
%

% some pipelines get super fancy, but why get super fancy when you have an iir filter, because
% moar infinite


nans=isnan(DATA);
DATA(nans)=0;

[b,a]=ellip(5,.2,40,[FC-BW/2 FC+BW/2]/(FS/2),'bandpass');

FILT_DATA=filtfilt(b,a,DATA);
FILT_DATA(nans)=nan;
