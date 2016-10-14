%%% slides for Max Planck
%

% load in xcorr data

%%

mat=kinect_read_csv([],',');


%%

[nrois,nsamples]=size(mat);
ts=1:nsamples;
ts=ts-299;
ts=ts/30;

%%

sortval=mean(mat(:,299:329),2);
[~,idx]=sort(sortval,'descend');


%%

% plot sorted matrix, clean up

figure();imagesc(ts,[],mat(idx,:));
ylimits=ylim();
hold on;
plot([0 0],[ylimits],'w-','linewidth',1.5);
plot([1 1],[ylimits],'w-','linewidth',1.5);
xlim([-5 5]);
set(gca,'TickDir','out','XTick',[-5:2.5:5],'TickLength',[0 0],'YTick',ylimits,'YTickLabel',[1 nrois],'FontSize',15);
ylabel('ROI');
xlabel('Time (s)');


%% 

% take top and bottom quartiles, plot average

use_rois=round(nrois*.25);

figure();
h(1)=plot(ts,mean(mat(idx(1:use_rois),:)),'r-')
hold on;
h(2)=plot(ts,mean(mat(idx(end-use_rois:end),:)),'b-');
ylim([-2 2]);
ylimits=ylim();
plot([0 0],ylimits,'w-');
xlim([-5 5]);
set(gca,'TickLength',[0 0],'XTick',[-5:2.5:5],'FontSize',15,'YTick',[-2 2]);
ylabel('Corr. (Z)');
xlabel('Time (s)');
L=legend(h,{'Top quartile','Bottom quartile'});
legend boxoff;
set(L,'FontSize',12);
