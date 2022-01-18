% t = Tiff(fullfile([folder_name_wr,'_rig.tif']),'r');
% M_rg = read(t);
% M_rg=load(fullfile([folder_name_wr,'_rig.tif']));
im=mean((M_rg),3);
im=im2double(im);
% im=normalize(im);
% im=(im-mean(im(:)))./std(im(:));
im=(im-min(im(:)))/(max(im(:))-min(im(:)));
im=imadjust(im,[0 0.2],[]);
% im=(im-min(im(:)))/(max(im(:))-min(im(:)));
figure;imshow(im);colormap(gray)
imwrite(im,fullfile([folder_name_wr,'average.tif']));

ui=input('Please press y to select a part of image: ','s');
if ui=='y'
    h=imfreehand;
    bw=createMask(h);
    im_n=bw.*im;
    im_n=imadjust(im_n,[min(im_n(im_n>0)),0.5],[]);
    figure;imshow(im_n);colormap(gray)
    imwrite(im_n,fullfile([folder_name_wr,'average_selected.tif']));
end

% im_max=max(im2double(M_rg),[],3);
% im_max=(im_max-min(im_max(:)))/(max(im_max(:))-min(im_max(:)));
% figure;imshow(im_max);colormap(gray)
% imwrite(im_max,fullfile([folder_name_wr,'max_proj.tif']));

% for i=1:10
%     im=mean(im2double(M_rg(:,:,501:1000)),3);
%     im=imadjust(im);
% im=(im-min(im(:)))/(max(im(:))-min(im(:)));
% figure;imshow(im);colormap(gray)

% for i=1:500%200:230;%1:500
%     imagesc(M_rg(:,:,i));title(num2str(i));colormap('gray');
% % imshow(imadjust(M_rg(:,:,i)));%colormap('gray');
% w = waitforbuttonpress;
% end
% 
% for i=700:730;%1:500
%     imagesc(M_rg(:,:,i));title(num2str(i));colormap('gray');
% % imshow(imadjust(M_rg(:,:,i)));%colormap('gray');
% w = waitforbuttonpress;
% end