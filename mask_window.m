
% Generation of Mask images 
%______________________________________________________________________________________________________________________________________________

function [mask1, mask2, mask3,mask4] = mask_window(win_size)

wsize = 9; 

rho=5;
theta=0;    % orizzontale
sigma2=5;
th_kernel=1e-3;
D=[rho, 0; 0, 1/rho];
R=[cos(theta), -sin(theta); sin(theta), cos(theta)];
S= inv( sigma2*R*D*R' );

[x,y]=meshgrid(-wsize:wsize, -wsize:wsize);
g1=zeros(size(x));

for k=1:length(x(:))
    g1(k)=(2*pi*sigma2)^-1 * exp(-0.5* [x(k), y(k)]*S*[x(k), y(k)]');
end

g2=zeros(size(x));
theta=0.5*pi;    % verticale
R=[cos(theta), -sin(theta); sin(theta), cos(theta)];
S= inv( sigma2*R*D*R' );
for k=1:length(x(:))
    g2(k)=(2*pi*sigma2)^-1 * exp(-0.5* [x(k), y(k)]*S*[x(k), y(k)]');
end

g3=zeros(size(x));
theta=0.25*pi;    % 45 gradi
R=[cos(theta), -sin(theta); sin(theta), cos(theta)];
S= inv( sigma2*R*D*R' );
for k=1:length(x(:))
    g3(k)=(2*pi*sigma2)^-1 * exp(-0.5* [x(k), y(k)]*S*[x(k), y(k)]');
end

g4=zeros(size(x));
theta=0.75*pi;    % 135 gradi
R=[cos(theta), -sin(theta); sin(theta), cos(theta)];
S= inv( sigma2*R*D*R' );
for k=1:length(x(:))
    g4(k)=(2*pi*sigma2)^-1 * exp(-0.5* [x(k), y(k)]*S*[x(k), y(k)]');
end


[x2,y2]=meshgrid(linspace(-wsize,wsize, 2*win_size+1), linspace(-wsize,wsize, 2*win_size+1));
g1=interp2(x, y, g1, x2, y2);
g2=interp2(x, y, g2, x2, y2);
g3=interp2(x, y, g3, x2, y2);
g4=interp2(x, y, g4, x2, y2);



mask1=g1>th_kernel;
mask2=g2>th_kernel;
mask3=g3>th_kernel;
mask4=g4>th_kernel;

% mask(:,:,1) = mask1; 
% mask(:,:,2) = mask2; 
% mask(:,:,3) = mask3; 
% mask(:,:,4) = mask4; 


end
%______________________________________________________________________________________________________________________________________________
