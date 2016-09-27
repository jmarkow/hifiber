%%%% script for visualizing the scalars-based analysis


% assume we've loaded in analysis_scalars.mat

scalar_names=fieldnames(obs_z);
camera_fs=30;
z_image=cellfun(@(x) obs_z.(x), scalar_names,'uni',0);
z_image=cat(2,z_image{:});
z_image_min=prctile(z_image(:),.5);
z_image_max=prctile(z_image(:),99.5);

% map z_image to indexed image (64 )

z_image=(z_image-z_image_min)./(z_image_max-z_image_min);
z_image=uint8(z_image*256);
z_image_rgb=ind2rgb(z_image',winter(256));

p_image_left=cellfun(@(x) obs_p.left.(x), scalar_names,'uni',0);
p_image_left=cat(2,p_image_left{:})';

p_image_right=cellfun(@(x) obs_p.right.(x), scalar_names,'uni',0);
p_image_right=cat(2,p_image_right{:})';

highlight=~(p_image_left<.01|p_image_right<.01);
[r,c]=find(highlight);

for i=1:length(r(:))
    z_image_rgb(r(i),c(i),:)=z_image_rgb(r(i),c(i),:).*.5;
end

time_vec=floor(size(z_image_rgb,2)/2);
time_vec=[-time_vec:time_vec]./camera_fs;

figure();image(time_vec,[],z_image_rgb);
set(gca,'YTick',[1:size(z_image_rgb,1)],'YTickLabels',regexprep(scalar_names,'_','\\_'));

% set(gca,'YTick',[1:length(scalar_names)],'YTickLabels',regexprep(scalar_names,'_','\\_'))


%figure();imagesc(cat(2,tmp{:})');
%set(gca,'YTick',[1:length(scalar_names)],'YTickLabels',regexprep(scalar_names,'_','\\_'))

% convert to rgb image, scale brightness by the p-value

%
% tmp=cellfun(@(x) obs_p.left.(x), scalar_names,'uni',0);
% figure();imagesc(-log10(cat(2,tmp{:}))');
% set(gca,'YTick',[1:length(scalar_names)],'YTickLabels',regexprep(scalar_names,'_','\\_'))
% caxis([-log10(.05) 3]);
%
% tmp=cellfun(@(x) obs_p.right.(x), scalar_names,'uni',0);
% figure();imagesc(-log10(cat(2,tmp{:}))');
% set(gca,'YTick',[1:length(scalar_names)],'YTickLabels',regexprep(scalar_names,'_','\\_'))
% caxis([-log10(.05) 3]);
