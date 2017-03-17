%%
%
%
%

%%

use_session=6;
use_mat_gcamp=extract_object(use_session).neural_data.analysis.syll_corr_gcamp;
use_mat_rcamp=extract_object(use_session).neural_data.analysis.syll_corr_rcamp;
use_mat_gcampz=zscore(use_mat_gcamp);
use_mat_rcampz=zscore(use_mat_rcamp);
use_mat_diff=use_mat_gcamp-use_mat_rcamp;

%%

%use_parameters=parameters{5};
ar_mat=params.ar_mat(:,:,1:30);
ar_mat=permute(ar_mat,[3 2 1]);
ar_mat=reshape(ar_mat,size(ar_mat,1)*size(ar_mat,2),[]);

[coef score]=pca(zscore(ar_mat)');
[coeff2 score2]=pca(use_mat_rcampz');
[coeff3 score3]=pca(use_mat_gcampz');

gcampmax=zeros(36,1);
for i=1:99
    tmp=corrcoef(score(idx(1:36),i),mean(use_mat_gcamp(60:90,:)));  
    gcampmax(i)=tmp(2,1);
end

rcampmax=zeros(36,1);
for i=1:99
    tmp=corrcoef(score(idx(1:36),i),mean(use_mat_rcamp(60:90,:)));  
    rcampmax(i)=tmp(2,1);
end

% corrmat=zeros(20,1);
% 
% for i=1:20
%     for j=1:20
%         tmp=corrcoef(score(idx(1:36),i),score2(:,j));
%         corrmat(i,j)=tmp(2,1);
%     end
% end
% 
% 
% corrmat2=zeros(20,20);
% 
% for i=1:20
%     for j=1:20
%         tmp=corrcoef(score(idx(1:36),i),score3(:,j));
%         corrmat2(i,j)=tmp(2,1);
%     end
% end
