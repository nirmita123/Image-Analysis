%% Single Image Haze Removal Using Dark Channel Prior

% I(x) = J(x)t(x) + A(1-t(x))
% I = observed intensity
% J = scene radiance                J(x)t(x)  = direct attenuation
%       Direct attenuation = scene radiance and delay in medium
% A = global atmospheric light      A(1-t(x)) = airlight
% t = medium transmission describing the portion of light that is not
% scattered and reaches the camera
%       Airlight = previously scattered light and shift of scene colors

clear all
clc

I = imread('city.jpg'); % Reading image
[xind, yind, z] = size(I); % Saving dimensions of image in xind, yind, z

figure;
imshow(I);

%% Generating Dark Channel Prior
% Finding minimum of RGB values of each pixel
for i = 1:xind
    for j = 1:yind
        Im(i,j)=min(I(i,j,:));
    end
end
% Applying min filter over patch size 15x15
ps = 15;
A = ordfilt2(Im,1,ones(ps,ps), 'symmetric');

figure;
imshow(A);

%% Haze Removal Using Dark Channel Prior
% Estimating the transmission
% t(x) = 1 - w min[over patch](min[over RGB] I/A)
% For this we would require A = atmospheric light
% Once we estimate the transmission, we will apply soft matting

%% Estimating the Atmospheric Light
% Finding top 0.1 percent brightest pixels in the dark channel
mat = reshape(A,[],1); % Reshape matrix into vector
B = sort(mat, 'descend'); % Sort the vector in descending order
B = B(1:0.001*xind*yind); % Changes B to top 0.1% values of B

%% Estimating the Transmission
% t(x) = 1 - w min[over patch](min[over RGB] I/A)
% Here we have to normalize each color channel independently in I/A
A1 = reshape(mat,xind,yind); % Reshape the matrix back to xind x yind
[row, col] = find(A1, B(1)); % Find the coordinates of the topmost value of B in A1
s = size(col); % Find the number of coordinates with maximum value
mazz = zeros(s(1), 3); % Create a 2D matrix sized s each row having 3 columns (RGB)
% We find image value corresponding to each row,col pair having maximum
% value in the airlight and store it in mazz
for i = 1:s(1)
    mazz(i, 1)=I(row(i), col(i), 1);
    mazz(i, 2)=I(row(i), col(i), 2);
    mazz(i, 3)=I(row(i), col(i), 3);
end
mazz = max(mazz); %  We find the row with maximum value in mazz
w = 0.95; % Since we do not wish to remove all haze from image

% Estimating the transmission
Ifa = [1]; % Declaring Ifa
% Converting to double for accurate calculation
Ifa = double(Ifa);
I = double(I);
mazz = double(mazz);
% Normalizing each color channel independently
for i = 1:xind
    for j = 1:yind
        Ifa(i,j)=min([I(i,j,1)/mazz(1), I(i,j,2)/mazz(2), I(i,j,3)/mazz(3)]);
    end
end
Im1=ordfilt2(Ifa,1,ones(ps,ps),'symmetric');
t = (1 - w.*Im1);

Impr = t*255;
Impr = uint8(Impr);

figure;
imshow(Impr);

%% Recovering the Scene Radiance
% Final scene radiance is J
% J(x) = (I(x)-A) / max(t(x),t0) + A
% Here t0 = 0.1
J = [1]; % Predefining J
for i = 1:xind
    for j = 1:yind
        J(i,j,1) = (I(i,j,1) - mazz(1))/max([t(i,j),0.1])+mazz(1);
        J(i,j,2) = (I(i,j,2) - mazz(2))/max([t(i,j),0.1])+mazz(2);
        J(i,j,3) = (I(i,j,3) - mazz(3))/max([t(i,j),0.1])+mazz(3);
    end
end
J = uint8(J); % Converting to uint8 for displaying image
figure;
imshow(J);
