clear;
clc;

%% read data from file 

filepath = './';
filename = 'CE_BSA_17';
ext      = '.tif';

% load txt image
data = importdata([filepath, filename, ext]);
data = double (max(data,0)); % clip negative values 
[nx,ny] = size(data);
N = 2^ceil(log2(max(nx,ny))+1); % size of image for Fourier transform, at least twice the size of a raw image

%% stich patches of the original image to form a periodic 'super-image' 
im = [data fliplr(data) data fliplr(data) data]; 
im = [flipud(im); im; flipud(im); im; flipud(im)]; 
%figure;imagesc(im);

% cut out the center part of the 'super-image' of size NxN 
[sx,sy] = size(im);
im = im(floor((sx-N)/2)+(1:N),floor((sy-N)/2)+(1:N)); 
%figure;imagesc(im);

% let the signal decay to zero towards the edges from R=N/2 to R=N
[X,Y] = meshgrid(1:N,1:N);
X = (X-N/2-1); Y = (Y-N/2-1); R = sqrt(X.^2+Y.^2); 
W = (2-4*R/N); W = min(W,1); W = max(W,0); 

im = W.*im;
%figure;imagesc(im);

%% Perform Fourier Reweighting
% Do the Fourier transform
fim = fftshift(fft2(im)); 
fa = abs(fim); % modulus
fp = angle(fim); % Phase angle

figure;imagesc(log(fa)); axis square; axis off

% calculate k-vectors of the Fourier-transformed image

fs = 1/0.045; % per µm, inverse of pixel size
k  = fs/2/(N/2).*R; 
km = 5.5; % km is 2NA / lambda_pump + 2NA / lambda_stokes
ind = double(k<km); % Step function to clip everything outside km to zero 
figure;imagesc(log(ind));axis square; axis off

%% coplot ft of the image and the frequency cutoff circle
figure;imagesc(log(fa)); hold on;

centers = [1024 1024];
radii = 400;
viscircles(centers,radii,'Color','r'); axis square; axis off; hold off

%% Define OTF and reweighting function
% 'Ideal' MTF
OTF = ind.*2/pi.*(acos(k./km)-k./km.*sqrt(1-(k./km).^2)); 
% Re-weighting function
ep = 0.25; 
w = ind./(OTF + ep.*k./km);
%figure;imagesc(w)
%%
% Back-transform the re-weighted FT 
tim = abs(ifft2(ifftshift(w.*fa.*exp(1i.*fp)))); 
% cut out the original field of the data 
tim = tim((N-nx)/2+(1:nx),(N-ny)/2+(1:ny)); 
tim = flipud(tim);
% display original data vs. re-weighted data 

figure;imagesc([data, tim])
colormap bone
%colormap(gray(256)) 
%caxis([0 350]) 
axis image 
axis off

%% Write image as tiff

opt_filepath = './';
ext = '.tif';
FR_out_filename = [filename, '_FR_git', ext];
imwrite(uint16(tim), [opt_filepath, FR_out_filename],'Compression','none');
%FR_out_filename = [filename, '_FR', ext];
%dlmwrite([opt_filepath, FR_out_filename], tim, 'delimiter','\t');
