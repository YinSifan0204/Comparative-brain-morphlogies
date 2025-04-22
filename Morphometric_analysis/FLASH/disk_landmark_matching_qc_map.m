function [map,s1,s2,f1_filled,f2_filled] = disk_landmark_matching_qc_map(v1,f1,lm1,v2,f2,lm2)

% This program computes a landmark-matching quasi-conformal parameterization
% of a simply-connected open brain surface onto the unit disk.
%  
% Both the initial sulcal landmark curve processing and the mapping 
% computation are included. 
%
% Input:
% v1: mx3 matrix for vertices of the source mesh (open surface)
% f1: px3 matrix for faces of the source mesh 
% lm1: kx1 cell array with each entry containing a row vector of indices of the landmark curves on the source mesh
% v2: nx3 matrix for vertices of the target mesh (open surface)
% f2: qx3 matrix for vertices of the target mesh
% lm2: kx1 cell array with each entry containing a row vector of indices of the landmark curves on the target mesh
%
% Output:
% map: mx2 coordinates of the disk landmark-matching quasi-conformal map
% s1: mx2 coordinates of the disk conformal parameterization of the source mesh
% s2: nx2 coordinates of the disk conformal parameterization of the target mesh
% f1_filled: px3 matrix for faces of the source mesh, 
%            possibly with holes filled to make it simply-connected
% f2_filled: qx3 matrix for faces of the target mesh, 
%            possibly with holes filled to make it simply-connected
%
% Written by Gary Pui-Tung Choi, 2023

if length(lm1)~=length(lm2) || isempty(lm1) || isempty(lm2)
    error('The number of landmark curves are different.');
end

if size(lm1,2) > size(lm1,1)
    lm1 = lm1';
end

if size(lm2,2) > size(lm2,1)
    lm2 = lm2';
end

bd1 = meshboundaries(f1);
bd2 = meshboundaries(f2);

if length(bd1) ~= 1
    warning(['Mesh 1 is not a simply-connected open surface! It contains ', ...
        num2str(length(bd1)), ...
        ' boundaries instead of 1. The extra holes will be filled.']);
    f1 = hole_filling(f1,1);
end

if length(bd2) ~= 1
    warning(['Mesh 2 is not a simply-connected open surface! It contains ', ...
        num2str(length(bd2)), ...
        ' boundaries instead of 1. The extra holes will be filled.']);
    f2 = hole_filling(f2,1);
end

f1_filled = f1;
f2_filled = f2;

% to avoid error in mobius transformation
if length(lm1) == 1
    lm1_temp = cell2mat(lm1);
    lm2_temp = cell2mat(lm2);
    
    lm1 = cell(2,1); 
    lm1{1} = lm1_temp(1:length(lm1_temp)/2);
    lm1{2} = lm1_temp(round((length(lm1_temp)+1)/2):end);
    
    lm2 = cell(2,1); 
    lm2{1} = lm2_temp(1:length(lm2_temp)/2);
    lm2{2} = lm2_temp(round((length(lm2_temp)+1)/2):end);
end

% ensure equal number of landmark points on each curve
[lm1a,lm2a,index] = process_lm(lm1,lm2);

fprintf('Landmark resampling completed.\n');

%% obtain the conformal parameterizations
fprintf('Computing the disk conformal map of the source mesh... ');
s1 = disk_conformal_map(v1,f1);
mu1 = beltrami_coefficient(s1,f1,v1);
mu1(abs(mu1)>0.9) = 0.9*exp(1i*angle(mu1(abs(mu1)>0.9)));
s1 = linear_beltrami_solver(s1,f1,mu1,bd1{1},s1(bd1{1},:));
[s1,~] =  mobius_area_correction_disk(v1,f1,s1);
fprintf('Done.\n');

fprintf('Computing the disk conformal map of the target mesh... ');
s2 = disk_conformal_map(v2,f2);
mu2 = beltrami_coefficient(s2,f2,v2);
mu2(abs(mu2)>0.9) = 0.9*exp(1i*angle(mu2(abs(mu2)>0.9)));
s2 = linear_beltrami_solver(s2,f2,mu2,bd2{1},s2(bd2{1},:));
[s2,~] =  mobius_area_correction_disk(v2,f2,s2);
fprintf('Done.\n');

fprintf('Initial landmark mismatch energy = %d\n', ...
    sum((s2(lm2a,1)-s1(lm1a,1)).^2+(s2(lm2a,2)-s1(lm1a,2)).^2));

%% Do an extra Mobius transformation to roughly align the landmarks

[s1,~] = mobius_landmark_alignment_disk(s1,lm1a,s2(lm2a,:));

fprintf('Landmark mismatch energy after Mobius alignment = %d\n', ...
    sum((s2(lm2a,1)-s1(lm1a,1)).^2+(s2(lm2a,2)-s1(lm1a,2)).^2));

%% Lamdmark-matching quasi-conformal map on the unit disk
% map = landmark_matching_qc_map_smooth_v3(s1,f1,lm1a',s2(lm2a',:),bd1{1},s1(bd1{1},:));
map = landmark_matching_qc_map_smooth_v4(v1,f1,s1,lm1a',s2(lm2a',:),bd1{1},s1(bd1{1},:));

fprintf('Final landmark mismatch energy = %d\n', ...
    sum((s2(lm2a,1)-map(lm1a,1)).^2+(s2(lm2a,2)-map(lm1a,2)).^2));
