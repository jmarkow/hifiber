%%
% assume we're looking at 15781 20160516131139

m=matfile('depth_masked.mat');
load('frame2_mask.mat');

% plot the second frame

use_data=m.depth_masked(:,:,1:100);
[r,c,nframes]=size(use_data);

[xx,yy]=meshgrid(1:c,1:r);

%%

tmp=double(use_data(:,:,2));
mask=mask&tmp>20;
use_features=[xx(mask) yy(mask) tmp(mask)];
all_features=[xx(:) yy(:) tmp(:)];

mu=mean(use_features);
sig=cov(use_features);

likelihood=mvnpdf(all_features,mu,sig);

tmp2=double(use_data(:,:,10));
likelihood2=mvnpdf([xx(:) yy(:) tmp2(:)],mu,sig);

%%

ax=[];
figure();
ax(1)=subplot(2,2,1);imagesc(tmp);
caxis([5 40]);
axis off;
ax(2)=subplot(2,2,2);imagesc(reshape(likelihood,size(tmp)));
caxis([1e-8 1e-4]);
axis off;
ax(3)=subplot(2,2,3);imagesc(tmp2);
caxis([5 40]);
axis off;
ax(4)=subplot(2,2,4);imagesc(reshape(likelihood2,size(tmp)));
caxis([1e-8 1e-4]);
linkaxes(ax,'xy');
axis off;
