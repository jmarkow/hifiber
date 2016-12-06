load('experiment_data_pca.mat','features');
load('pca_analysis.mat','lag_mat','rnd_summary');

zmu_green=zeros(length(lag_mat),size(lag_mat(1).proc1,2));
zmu_red=zeros(length(lag_mat),size(lag_mat(1).proc2,2));
% pmu_right=zeros(size(zmu));
% pmu_left=zeros(size(zmu));
mu_green=zeros(size(zmu_green));
mu_red=zeros(size(zmu_red));

for i=1:length(lag_mat)
   
   % express as X fold over the phase rnd
   
   mu_green(i,:)=mean(lag_mat(i).proc1);
   mu_red(i,:)=mean(lag_mat(i).proc2);
   
   zmu_green(i,:)=(mu_green(i,:)-mean(rnd_summary(i).proc1))./std(rnd_summary(i).proc1);
   zmu_red(i,:)=(mu_red(i,:)-mean(rnd_summary(i).proc2))./std(rnd_summary(i).proc2);
   % cut off by p-value?
   
   % exclude all significant pcs, two-tailed .05
   
   pmu_left_green(i,:)=mean(repmat(mu_green(i,:),[size(rnd_summary(i).proc1,1) 1])>rnd_summary(i).proc1);
   pmu_right_green(i,:)=mean(repmat(mu_green(i,:),[size(rnd_summary(i).proc1,1) 1])<rnd_summary(i).proc1);
      
end

% only include pcs with p-values<.01 in reasonable lags (-1/+1sec)
% 
% alpha_cut=.05;
% right_idx=find(sum(pmu_right(:,30:150)'<=alpha_cut)>0);
% left_idx=find(sum(pmu_left(:,30:150)'<=alpha_cut)>0);
% 
% zmu_right=zmu(right_idx,:);
% zmu_left=zmu(left_idx,:);