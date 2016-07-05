function [WINS WINS_T]=kinect_analysis_win_data(LOCS,WIN_SIZE,DATA)
%
%
%
%
%
%
%

%

WIN_T=[-WIN_SIZE:WIN_SIZE];
WINS=zeros(length(WIN_T),length(LOCS));
to_del=[];

for i=1:length(LOCS)

  left_edge=LOCS(i)-WIN_SIZE;
  right_edge=LOCS(i)+WIN_SIZE;

  if left_edge>0 & right_edge<size(DATA,1)
    WINS(:,i)=DATA(left_edge:right_edge);
  else
    to_del=[to_del i];
  end
end

WINS(:,to_del)=[];
