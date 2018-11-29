clear all
clc
% Reading image
I = imread('city.jpg');
[xind, yind, z] = size(I);
% Finding minimum of RGB values of each pixel
for i = 1:xind
    for j = 1:yind
        Im(i,j)=min(I(i,j,:));
    end
end
% Patch size
ps = 15;
% Minimum filter on patch size of 15x15
A = ordfilt2(Im,1,ones(ps,ps), 'symmetric');
% Estimating the atmospheric light
% Finding top 0.1 percent brightest pixels in the dark channel
mat = reshape(A,[],1);
B = sort(mat, 'descend');
B = B(1:0.001*xind*yind);
A1 = reshape(mat,xind,yind);
[row, col] = find(A1, B(1));
s = size(col);
mazz = zeros(s(1), 3);
for i = 1:s(1)
    mazz(i, 1)=I(row(i), col(i), 1);
    mazz(i, 2)=I(row(i), col(i), 2);
    mazz(i, 3)=I(row(i), col(i), 3);
end
mazz = max(mazz);
% w defined as 0.95 in the paper
w = 0.95;
% Estimating the transmission
Ifa = [1];
Ifa = double(Ifa);
I = double(I);
mazz = double(mazz);
for i = 1:xind
    for j = 1:yind
        Ifa(i,j)=min([I(i,j,1)/mazz(1), I(i,j,2)/mazz(2), I(i,j,3)/mazz(3)]);
    end
end
%Ifa = Ifa.*255;
Im1=ordfilt2(Ifa,1,ones(ps,ps),'symmetric');
t = (1 - w.*Im1);

% Converting to uint8 to display image
Ifa = uint8(Ifa);
%figure;
%imshow(t);
J = [1];
% Recovering the scene radiance
for i = 1:xind
    for j = 1:yind
        J(i,j,1) = (I(i,j,1) - mazz(1))/max([t(i,j),0.1])+mazz(1);
        J(i,j,2) = (I(i,j,2) - mazz(2))/max([t(i,j),0.1])+mazz(2);
        J(i,j,3) = (I(i,j,3) - mazz(3))/max([t(i,j),0.1])+mazz(3);
    end
end
I = uint8(I);
figure;
imshow(I);
J = uint8(J);
figure;
imshow(J);
