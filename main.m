% This script comuptes the invariants of two images on the sphere. The
% second image is being rotated in SO(3) using D-Wigner matrices.
% Tamir Bendory, Dec 2017

close all;

n = 10;  %must be even for chebfun (we need to find a faster way to compute SH!)

% If 1, the shift and in-plane rotation takes place in the image space. 
%Otherwise, it is perfomed after the projection on the SH coefficients on the sphere 
shift_space_domain = 0; 

if shift_space_domain
rot_angle = 0; % in-plane rotation in image space
shift = [1,0]; % shift in image space
zp = max(shift); % zero padding of the image
else 
    zp = 0;
end

% loading the example image -- "Cat's Eye Nebula"
rgb = imread('ngc6543a.jpg'); I = rgb2gray(rgb);
I = I(150:450,150:450);
indices = floor(linspace(1,size(I,1),n-zp));
I = double(imresize(I,[n-2*zp,n-2*zp]));
func = zeros(n); func(zp+1:n-zp,zp+1:n-zp) = I;
%func = func/(max(max(func)));

%figure; subplot(121); imshow(func); title('original image');
%subplot(122); imshow(func_rot); title('rotated image');

% moving to polar coordinates
theta_dom = linspace(-pi,pi,n);
h = (pi/4)/(n-1); radial_dom = 0:h:(pi/4);
phi = 0:h:2*pi;
[M, K] = meshgrid(theta_dom,radial_dom);
[xx, yy] = pol2cart(M, K);
x = linspace(-.5,.5,n);
y = linspace(-.5,.5,n);
[X,Y] = meshgrid(x,y);
[Az, El] = meshgrid(theta_dom, phi);

ext_func = zeros(size(Az));
ext_func(1:n,1:n) = func;
% projection onto the sphere
F = spherefun(ext_func);
SH_matrix = compute_SH_coeff(F,n);

if shift_space_domain
    func_rot = imrotate(func,rot_angle,'crop','bicubic');
    func_rot = circshift(func_rot,shift);
    ext_func_rot = zeros(size(Az));
    ext_func_rot(1:n,1:n) = func_rot;
    F_rot = spherefun(ext_func_rot);
    SH_matrix_rot = compute_SH_coeff(F_rot,n);
else
    rot_angle = 2*pi*rand(3,1);% random rotation in SO(3)
    SH_matrix_rot = cell(n+1,1);
    SH_matrix_rot{1}= SH_matrix{1};
    for j=1:n
        SH_matrix_rot{j+1} = wigner_d(j, rot_angle)*SH_matrix{j+1};
    end
end

%% computing invariants
[M1,M2,M3] = compute_invariants(SH_matrix);
[M1_rot,M2_rot,M3_rot] = compute_invariants(SH_matrix_rot);

%% computing error

err1 = abs(M1-M1_rot)/abs(M1);
err2 = norm(M2(:)-M2_rot(:))/norm(M2(:));
err3 = norm(M3(:)-M3_rot(:))/norm(M3(:));

display(strcat('err1 = ',num2str(err1)));
display(strcat('err2 = ',num2str(err2)));
display(strcat('err3 = ',num2str(err3)));
