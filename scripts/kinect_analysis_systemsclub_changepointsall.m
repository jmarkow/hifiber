%
%
%
%

files={};
files{1}=('/Users/jmarkow/Desktop/workspace/photometry/15781/analysis/delta_score.mat');
files{2}=('/Users/jmarkow/Desktop/workspace/photometry/15783/analysis/delta_score.mat');
files{3}=('/Users/jmarkow/Desktop/workspace/photometry/15784/analysis/delta_score.mat');

delta_score_all={};

for i=1:length(files)
    load(files{1},'delta_score')
    
    for j=1:length(delta_score)
        delta_score_all{end+1}=delta_score{j};
    end
    
end

%%


ibi=[];

for i=1:length(delta_score_all)
    
    tmp=delta_score_all{i};
    thresh=mean(tmp)+.1*std(tmp);
    [~,tmp2]=findpeaks(tmp,'minpeakheight',thresh);
    ibi=[ibi;diff(tmp2)'];
    
end