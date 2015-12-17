function [d] = co_hog_des(img)

[m, n, k] = size(img);
%img = img(floor(2*m/3):m,:,:);

img = imResample(single(img),[256 256])/255;
H = hog(img, 64, 16);
d = reshape(H, [1024,1,1]);