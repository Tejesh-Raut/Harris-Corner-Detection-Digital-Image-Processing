function myHarrisCornerDetector( sigma1, sigma2, k )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Referred algorithm from http://www.cse.psu.edu/~rtc12/CSE486/lecture06.pdf
inpImg = '..\data\boat.mat';
f1 = load(inpImg,'-mat');
f2 = f1.imageOrig;
original = mat2gray(f2);
min1 = min(min(original));
max1 = max(max(original));
original = ((original - min1)*1/(max1 - min1));

gfilter = fspecial('gaussian', 5, sigma1);
[dxmask, dymask] = gradient(gfilter);
%dxmask = [-1, 0, 1; -2, 0, 2; -1, 0, 1];
%dymask = [1, 2, 1; 0, 0, 0; -1, -2, -1];

derivativeX = conv2(original, dxmask, 'same');
derivativeY = conv2(original, dymask, 'same');

%derivativeX = ((derivativeX - min(min(derivativeX)))*1/(max(max(derivativeX)) - min(min(derivativeX))));
%derivativeY = ((derivativeY - min(min(derivativeY)))*1/(max(max(derivativeY)) - min(min(derivativeY))));

derivativeXX = derivativeX.*derivativeX;
derivativeYY = derivativeY.*derivativeY;
derivativeXY = derivativeX.*derivativeY;

g = fspecial('gaussian', 5, sigma2);
Sx2 = conv2(derivativeXX, g, 'same');
Sy2 = conv2(derivativeYY, g, 'same');
Sxy = conv2(derivativeXY, g, 'same');

[Srows, Scols ]= size(Sx2);
R = zeros(Srows, Scols);
harrisCornerness = R;
eigenvalue1 = zeros(Srows, Scols);
eigenvalue2 = zeros(Srows, Scols);
for i = 1:Srows
    for j = 1:Scols
        H = [Sx2(i,j), Sxy(i,j); Sxy(i,j), Sy2(i,j)];
        R(i, j)  = det(H) - (k*trace(H)*trace(H));
        harrisCornerness(i, j) = R(i, j);
        eigenvec = eig(H);
        eigenvalue1(i, j) = eigenvec(1);
        eigenvalue2(i, j) = eigenvec(2);
    end;
end;
%disp(eigenvalue1);
%disp(eigenvalue2);
myNumOfColors = 200;
myColorScale = [[0:1/(myNumOfColors - 1):1]',[0:1/(myNumOfColors - 1):1]' , [0:1/(myNumOfColors - 1):1]' ];

disp(' Tuned parameter values are: ');
disp([' Sigma1 is ',num2str(sigma1)]);
disp([' Sigma2 is ',num2str(sigma2)]);
disp([' K is ',num2str(k)]);

figure(1);
subplot(1, 2, 1);
imagesc ((derivativeX));
title('Derivative in X direction');
colormap (myColorScale);
colormap (jet);
daspect ([1 1 1]);
axis tight;
colorbar
subplot(1, 2, 2);
imagesc ((derivativeY));
title('Derivative in Y direction');
colormap (myColorScale);
colormap (jet);
daspect ([1 1 1]);
axis tight;
colorbar

figure(2);
subplot(1, 2, 1);
imagesc ((eigenvalue1));
title('First eigenvalue of the structure tensor');
colormap (myColorScale);
colormap (jet);
daspect ([1 1 1]);
axis tight;
colorbar
subplot(1, 2, 2);
imagesc ((eigenvalue2));
title('Second eigenvalue of the structure tensor');
colormap (myColorScale);
colormap (jet);
daspect ([1 1 1]);
axis tight;
colorbar

figure(3);
subplot(1, 1, 1);
imagesc ((harrisCornerness));
title('Harris corner-ness measure');
colormap (myColorScale);
colormap (jet);
daspect ([1 1 1]);
axis tight;
colorbar

imwrite(original,'..\images\boat.png');
imwrite(derivativeX,'..\images\boatDerivativeX.png');
imwrite(derivativeY,'..\images\boatDerivativeY.png');
imwrite(eigenvalue1,'..\images\boatEigenvalue1.png');
imwrite(eigenvalue2,'..\images\boatEigenvalue2.png');
imwrite(harrisCornerness,'..\images\boatHarrisCornerness.png');

end

