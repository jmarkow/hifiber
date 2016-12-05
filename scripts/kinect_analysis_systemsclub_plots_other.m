%% all changepoints

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

ibi=[];

for i=1:length(delta_score_all)
    
    tmp=delta_score_all{i};
    thresh=mean(tmp)+.1*std(tmp);
    [~,tmp2]=findpeaks(tmp,'minpeakheight',thresh);
    ibi=[ibi;diff(tmp2)'];
    
end

%%

camera_fs=30;

y=histc(ibi*1/camera_fs,bins);
ci=bootci(1e3,{@(x) histc(x*1/camera_fs,bins),ibi},'type','per','alpha',.01);

%%

figs.hist_fig=figure('PaperPositionMode','auto','Position',[100 100 400 300]);
whitebg(figs.hist_fig);
bins=[0:.05:1.4]*1e3;
g=gramm('x',(ibi*1/camera_fs)*1e3);
g.stat_bin('geom','stairs','edges',bins,'normalization','probability');
g.draw();
set(g.facet_axes_handles,'TickLength',[.015 .015],'TickDir','out','YLim',[-.01 .25],'xlim',[0 1400])
offsetAxes(g.facet_axes_handles)
set(figs.hist_fig,'color',[0 0 0],'InvertHardcopy','off')
markolab_multi_fig_save(figs.hist_fig,'~/Desktop/quickfigs','systemclub_changepoint_hist','eps,png,fig','renderer','painters');
