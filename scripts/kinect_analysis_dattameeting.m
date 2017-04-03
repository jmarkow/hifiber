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
ar_mat=params.ar_mat(beh(1).states(idx)+1,:,end-20:end-1);
ar_mat=permute(ar_mat,[3 2 1]);
ar_mat=reshape(ar_mat,size(ar_mat,1)*size(ar_mat,2),[]);

[coef score]=pca(zscore(ar_mat)');
[coeff2 score2]=pca(use_mat_rcampz');
[coeff3 score3]=pca(use_mat_gcampz');


%%

rsz=max(use_mat_rcamp(60:150,:));
rsz_max=prctile(rsz,95);
rsz_min=prctile(rsz,5);
rsz=(rsz-rsz_min)./(rsz_max-rsz_min);
rsz(rsz<0)=0;
rsz(rsz>1)=1;

gsz=max(use_mat_gcamp(60:150,:));
gsz_max=prctile(gsz,80);
gsz_min=prctile(gsz,1);
gsz=(gsz-gsz_min)./(gsz_max-gsz_min);
gsz(gsz<0)=0;
gsz(gsz>1)=1;

% gsz=max(use_mat_gcamp(60:150,useidx));
% gsz=(gsz-(gsz))./(max(gsz)-min(gsz));
gsz=min(gsz*3,1);
rsz=min(rsz*3,1);

% gsz=max(gsz,.2);
% rsz=max(rsz,.2);

%%



%mapped_data2=tsne(ar_mat(:,idx(1:40))',[],3,50,30);
figure(); scatter3(mapped_data2(:,1),mapped_data2(:,2),mapped_data2(:,3),100,[gsz(1:50)' rsz(1:50)' zeros(size(rsz(1:50)')) ],'filled')
%%

%%

num=5;
examples=10;
list=[32 83 11 93 16 37];
%list=24;
examples=10;
use_frames=extract_object(num).load_oriented_frames(true);

for i=1:length(list)
    
    [matches,score]=extract_object(num).template_match(list(i),1,30);
    pk_vals=score(matches);
    [~,pk_sorting]=sort(pk_vals,'descend');
    matches_sorted=matches(pk_sorting);
    
    for j=1:min(examples,length(matches_sorted))
    
        hit_onset=matches_sorted(j);
        syll_idx=beh(num).labels==list(i);
        syll_idx=score(hit_onset-60:hit_onset+60)>1;
        mov_frames=use_frames(:,:,hit_onset-60:hit_onset+60);
        marker_coords=cell(1,size(mov_frames,3));
        
        for k=1:length(marker_coords)
            if syll_idx(k)
                marker_coords{k}=[(1:10)' (1:10)'];
            end
        end
        
        kinect_extract.animate_direct(mov_frames,...
            'cmap',jet(256),'marker_coords',marker_coords,'filename',sprintf('syllable_%i_%i',list(i),j),...
            'clim',[0 80]);
        
    end
    
end

%%


gcampmax=zeros(95,1);
for i=1:95
    tmp=corrcoef(score(idx(1:size(use_mat_gcamp,2)),i),mean(use_mat_gcamp(60:90,:)));  
    gcampmax(i)=tmp(2,1);
end

rcampmax=zeros(95,1);
for i=1:95
    tmp=corrcoef(score(idx(1:size(use_mat_gcamp,2)),i),mean(use_mat_rcamp(60:90,:)));  
    rcampmax(i)=tmp(2,1);
end

%%

% goddamn I'm scaling this shit



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
