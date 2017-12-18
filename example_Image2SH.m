% script name: "example_Image2SH"
%
% In this script we demonstrate on a simple example image the process of 
% taking an image (ideally given in polar coordnates, here we actually 
% interpolate it to polar coordinates) and express it in fourier 
% coodinates on the sphere (Spherical Harminics or SH) after projecting it.
% 
% For simplicity, the SH coordinates are calculated using Chebfun ("http://www.chebfun.org") according to
% "Computing with functions in spherical and polar geometries I. The sphere"
%
% NS, Dec 2017

% In this example we take a standard image and move it to polar coordinates

n = 40;  %must be even for chebfun
close all;
tic

% loading the example image -- "Cat's Eye Nebula"
rgb = imread('ngc6543a.jpg');
I = rgb2gray(rgb);
figure; imshow(I); title('Original image');
I = I(150:450,150:450);
indices = floor(linspace(1,size(I,1),n));
func = double(I(indices,indices));
func = func/(max(max(func)));

% moving to polar coordinates
theta_dom = linspace(-pi,pi,n);
h = (pi/4)/(n-1); radial_dom = 0:h:(pi/4);
phi = 0:h:2*pi;
[M, K] = meshgrid(theta_dom,radial_dom);
[xx, yy] = pol2cart(M, K);
x = linspace(-.5,.5,n);
y = linspace(-.5,.5,n);
[X,Y] = meshgrid(x,y);

% interpolate pixel values for unknwon coordinates and present
% basically, this is the place to do 
%   it all over again, e.g. eith "Guassian grid"
out = interp2(X, Y, single(func), xx, yy);
figure; surf(X, Y ,func); view(2); 
figure; surf(xx, yy ,out); view(2); 

% mapping and its illustration 
[Az, El] = meshgrid(theta_dom, phi);

ext_func = zeros(size(Az));
ext_func(1:n,1:n) = func;

D_x = cos(Az).*sin(El); 
D_y = sin(Az).*sin(El);
D_z = cos(El); 

figure;
Hm = surf(D_x, D_y, D_z, ext_func);
set(Hm, 'EdgeAlpha', 0.1)
title('The discretization on the sphere');

% define the function on the equidistant grid
F = spherefun(ext_func);
figure;
plot(F), title('The image on the sphere'), colorbar, axis off

% calculating the coordinates
degree = 40; k = 1;
coeffs = zeros((degree+1)^2,3);
for l = 0:degree
    for m = -l:l
        Y = spherefun.sphharm(l,m);
        coeffs(k,1) = sum2(F.*Y);
        coeffs(k,2:3) = [l m];
        k = k + 1;
    end
end

figure;
st = stem3(coeffs(:,2),coeffs(:,3),abs(coeffs(:,1)),'filled'); 
set(st,'LineWidth',.4)
ylim([-degree degree])
set(gca,'ZScale','log'), set(gca,'Xdir','reverse'), view([-13 18])
xlabel('$\ell$','Interpreter','Latex'), ylabel('m'), zlabel('|coeffs|')
set(gcf, 'Position', get(0,'ScreenSize'));

% the projection
fproj = spherefun([]);
k = 1;
for l = 0:(degree)
    for m = -l:l
        fproj = fproj + coeffs(k,1)*spherefun.sphharm(l,m);
        k = k + 1;
    end
end

figure;
plot(fproj), title(['Degree ',num2str(degree),' spherical harmonic projection'])
colorbar, axis off

figure;
plot(F-fproj), title('Error in the spherical harmonic projection')
colorbar, colorbar, axis off
fprintf('The relative error is %f\n ', norm(F-fproj)/norm(F));
fprintf('Runtime is %d\n ', floor(toc));

