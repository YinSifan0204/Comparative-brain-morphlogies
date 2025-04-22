function [map_mobius_disk,x] =  mobius_landmark_alignment_disk(map,landmark,target)

% Find an optimal Mobius transformation for reducing the landmark mistmatch
% of a disk conformal parameterization.
%
% Input:
% map: nv x 2 vertex coordinates of the disk conformal parameterization
% landmark: k x 1 vertex indices of the landmarks
% target: k x 2 coordinates of the target positions
% 
% Output:
% map_mobius_disk: nv x 2 vertex coordinates of the updated disk conformal parameterization
% x: the optimal parameters for the Mobius transformation, where
%    f(z) = e^(it)*\frac{z-a}{1-\bar{a} z}
%    x(1): |a| (0 ~ 1)        magnitude of a
%    x(2): arg(a) (0 ~ 2pi)   argument of a
%    x(3): t (0 ~ 2pi)   rotation
%
% Written by Gary Pui-Tung Choi, 2023

%%
z = complex(map(:,1),map(:,2));

if size(target,2) == 2
    w = complex(target(:,1),target(:,2));
else
    w = target;
end

% objective function: difference in landmark position under the Mobius
% transformation
d_mismatch = @(x) sum(abs(...
    exp(1i*x(3))*(z(landmark)-x(1)*exp(1i*x(2)))./(1-conj(x(1)*exp(1i*x(2)))*z(landmark))...
    - w).^2);


% Optimization setup
x0 = [0,0,0]; % initial guess, try something diferent if the result is not good
lb = [0,-2*pi,-2*pi]; % lower bound for the parameters
ub = [1,2*pi,2*pi]; % upper bound for the parameters
options = optimoptions('fmincon','Display','off','MaxIter',150);

% Optimization (may further supply gradients for better result, not yet implemented)
x = fmincon(d_mismatch,x0,[],[],[],[],lb,ub,[],options);

% obtain the conformal parameterization with area distortion corrected
fz = exp(1i*x(3))*(z-x(1)*exp(1i*x(2)))./(1-conj(x(1)*exp(1i*x(2)))*z);
map_mobius_disk = [real(fz), imag(fz)];

end
