function I3=prepro(image)

imgBlur = image;%imread('riceblurred.png');
sharpCoeff = [0 0 0;0 1 0;0 0 0]-fspecial('laplacian',1);

imgSharp = imfilter(imgBlur,sharpCoeff,'symmetric');
% figure
% imshowpair(image,imgSharp,'montage')
% title('Original Image and Sharpened Image')

I=image;%imgSharp;
se = strel('disk',10);  %70
background=imerode(I,se);
I2=imreconstruct(background,I);
I3=I-imreconstruct(background,I);
% figure;
% subplot(2,2,1);imagesc(I);
% title('Original Etna Image')
% pbaspect([2 2 2])
% colorbar
% subplot(2,2,2);imagesc(background);
% title('Background Image')
% pbaspect([1 1 1])
% colorbar
% subplot(2,2,3);imagesc(I2);
% title('Reconstructed Image')
% pbaspect([1 1 1])
% colorbar
% subplot(2,2,4);imagesc(I3);
% title('TAB Image')
% pbaspect([1 1 1])
% colorbar
% colormap gray
