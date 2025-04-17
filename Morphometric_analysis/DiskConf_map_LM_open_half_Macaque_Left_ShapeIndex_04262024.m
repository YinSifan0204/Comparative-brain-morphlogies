% Compare real and simulated brains using disk conformal mapping parameterization
%           
%                              Sifan Yin 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear all; close all;set(0,'defaultfigurecolor','w')

addpath(genpath('Macaque_adult_open_half'));
parentDirectory = fileparts(pwd);
mpath = [parentDirectory,'/FLASH'];
addpath(mpath)

i = 1;
lr = {'left','right'};

name = ['Macaque_adult_',lr{i}];
new_mat = [name,'.mat'];

load(new_mat);

v3 = vertex_real;
v1 = vertex_simu;
v2 = vertex_gel;

f3 = face_real;
f1 = face_simu;
f2 = face_gel;


curvature1 = tricurv(f1,v1);
mean_curv1 = curvature1.km;
Gauss_curv1 = curvature1.kg;
ShapeIndexHK1 = 2/pi*atan(mean_curv1./sqrt(mean_curv1.^2-Gauss_curv1));

curvature2 = tricurv(f2,v2);
mean_curv2 = curvature2.km;
Gauss_curv2 = curvature2.kg;
ShapeIndexHK2 = 2/pi*atan(mean_curv2./sqrt(mean_curv2.^2-Gauss_curv2));

curvature3= tricurv(f3,v3);
mean_curv3 = curvature3.km;
Gauss_curv3 = curvature3.kg;
ShapeIndexHK3 = 2/pi*atan(mean_curv3./sqrt(mean_curv3.^2-Gauss_curv3));

index = find(abs(mean_curv2)>100);
mean_curv2(index) = 0;
index = find(isnan(mean_curv2));
mean_curv2(index) = 0;

mean(mean_curv2)
std(mean_curv2)

%%%%%%%%%%%%%%%%%%%
figure(10)
binWidth = 0.01;
Edges = -1:binWidth:1+binWidth;
num_vertex1  = length(ShapeIndexHK1);
num_vertex2  = length(ShapeIndexHK2);
num_vertex3  = length(ShapeIndexHK3);
centers = Edges(1:end-1)+binWidth/2;

counts1 = histcounts(ShapeIndexHK1,Edges)/num_vertex1;
counts2 = histcounts(ShapeIndexHK2,Edges)/num_vertex2;
counts3 = histcounts(ShapeIndexHK3,Edges)/num_vertex3;

plot(centers,counts3,'k.',centers,counts1,'b.',centers,counts2,'r.','MarkerSize',15);
xlim([-1,1])
xlabel('Shape Index')
ylabel('Probability')
title('Macaque')
%legend('real','simulated','gel','Location','northwest');
set(gca,'Fontsize',20,'Fontname','Arial')

set(0,'defaultlineLinewidth',3);
set(0,'defaultaxesLinewidth',3);
set(0,'defaultaxesfontsize',22);
print('Compare Three shape index_macaque_left_scatter.tif','-dtiff')

%%
figure(11)
I1 = threshold(mean_curv1,1);
I2 = threshold(mean_curv2,1);
I3 = threshold(mean_curv3,1);
% rescale both curvatures to [0,1]
I1 = (I1-min(I1))/(max(I1)-min(I1));
I2 = (I2-min(I2))/(max(I2)-min(I2));
I3 = (I3-min(I3))/(max(I3)-min(I3));

binWidth = 0.005
Edges = 0:binWidth:1+binWidth;
num_vertex1  = length(ShapeIndexHK1);
num_vertex2  = length(ShapeIndexHK2);
num_vertex3  = length(ShapeIndexHK3);
centers = Edges(1:end-1)+binWidth/2;

counts1 = histcounts(I1,Edges)/num_vertex1;
counts2 = histcounts(I2,Edges)/num_vertex2;
counts3 = histcounts(I3,Edges)/num_vertex3;

plot(centers,counts3,'k.',centers,counts1,'b.',centers,counts2,'r.');
xlim([0,1])
ylim([0,0.02])
xlabel('Rescaled mean curvature')
ylabel('Probability')
title('Mean curvature distribution')
legend('real','simulated','gel','Location','northwest');
set(gca,'Fontsize',18,'Fontname','Times New Roman')
set(0,'defaultlineLinewidth',1.5);
set(0,'defaultaxesLinewidth',1.5);
set(0,'defaultaxesfontsize',22);
print('Compare Three meancurv_macaque_left_scatter.tif','-dtiff')

%% Real ferret jet plot
figure(12)
set(gcf,'position',[50,50,1000,800])
%plot_meshg_lm(v3,f3,lm3,ShapeIndex3);
plot_meshj(v3,f3,ShapeIndexHK3);
%title([lr{i},' real ferret brain'])
view([0 10]); caxis([-1,1]) 
c = colorbar;set(c,'YTick',[-1:0.5:1])
print([name,'_real_lm_SI_jet_left.tif'],'-dtiff');


%% With landmarks
landmark_real = [name,'_real_lm'];
load(landmark_real)
lm_real = LM';
lm3 = lm_real;

landmark_simu = [name,'_simu_lm']
load(landmark_simu)
lm_simu = LM';
lm1 = lm_simu;

landmark_gel = [name,'_gel_lm']
load(landmark_gel)
lm_gel = LM';
lm2 = lm_gel;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1) 
set(gcf,'position',[50,50,1000,800])
plot_meshg_lm(v3,f3,lm3,threshold(mean_curv3,0.75));
%title([lr{i},' real ferret brain'])
view([0 10])  
print([name,'_curv_real_lm.tif'],'-dtiff');

figure(101) 
set(gcf,'position',[50,50,1000,800])
plot_meshg_lm(v3,f3,lm3,ShapeIndexHK3);
%title([lr{i},' real ferret brain'])
view([0 10])  
caxis([-1,1])
print([name,'_ShapeIndex_real_lm.tif'],'-dtiff');




% Rotate simulated Macaque
v1xz = v1(:,[1,3]);
theta = 10;
theta = theta/180*pi;
RotMat = [cos(theta), -sin(theta);sin(theta),cos(theta)];

v1xz_rot = v1xz*RotMat';
v1_rot = zeros(size(v1));
v1_rot(:,2) = v1(:,2);
v1_rot(:,[1,3]) = v1xz_rot;
v1 = v1_rot;


figure(2)
set(gcf,'position',[50,50,1000,800])
plot_meshb_lm(v1,f1,lm1,threshold(mean_curv1,0.75));
%title([lr{i},' simulated ferret brain'])
view([0 10])  
print([name,'_curv_simu_lm.tif'],'-dtiff');

figure(102)
set(gcf,'position',[50,50,1000,800])
plot_meshb_lm(v1,f1,lm1,ShapeIndexHK1);
%title([lr{i},' simulated ferret brain'])
view([0 10])  
caxis([-1,1])
print([name,'_ShapeIndex_simu_lm.tif'],'-dtiff')




v2xz = v2(:,[1,3]);
theta = -10;
theta = theta/180*pi;
RotMat = [cos(theta), -sin(theta);sin(theta),cos(theta)];
v2xz_rot = v2xz*RotMat';
v2_rot = zeros(size(v2));
v2_rot(:,2) = v2(:,2);
v2_rot(:,[1,3]) = v2xz_rot;
v2 = v2_rot;

figure(3) 
set(gcf,'position',[50,50,1000,800])
plot_meshp_lm(v2,f2,lm2,threshold(mean_curv2 ,0.3));
%title([lr{i},' real ferret brain'])
view([0 10])  
print([name,'_curv_gel_lm.tif'],'-dtiff');
figure(103) 
set(gcf,'position',[50,50,1000,800])
plot_meshp_lm(v2,f2,lm2,ShapeIndexHK2);
%title([lr{i},' real ferret brain'])
view([0 10])   
caxis([-1,1])
print([name,'_ShapeIndex_gel_lm.tif'],'-dtiff');










%% Disk map for a triple of surfaces

% map13: disk parameterization of Mesh 1 based on the landmark
%        positions of Mesh 3
%
% map23: disk parameterization of Mesh 2 based on the landmark
%        positions of Mesh 3

[map13,map23,s1,s2,s3,f1_filled,f2_filled,f3_filled] = ...
    disk_landmark_matching_qc_map_triple1029(v1,f1,lm1,v2,f2,lm2,v3,f3,lm3);


% [map13,map23,s1,s2,s3,f1_filled,f2_filled,f3_filled] = ...
%     disk_landmark_matching_qc_map_triple(v1,f1,lm1,v2,f2,lm2,v3,f3,lm3);

figure(4) 
set(gcf,'position',[50,50,1000,800])
plot_meshg_lm_disk(s3,f3_filled,lm3,threshold(mean_curv3,0.75));
view([90 90])  
print([name,'_mapping_real_lm.tif'],'-dtiff');
%title([lr{i},' real ferret brain'])
figure(41) 
set(gcf,'position',[50,50,1000,800])
plot_meshg_lm_disk(s3,f3_filled,lm3,ShapeIndexHK3);
view([90 90])  
print([name,'_mapping_real_lm_ShapeIndex.tif'],'-dtiff');
%title([lr{i},' real ferret brain'])

figure(5)
set(gcf,'position',[50,50,1000,800])
plot_meshb_lm_disk(map13,f1_filled,lm1,threshold(mean_curv1,0.75));
%title([lr{i},' simulated ferret brain'])
view([90 90])  
print([name,'_mapping_simu_lm.tif'],'-dtiff');
figure(51)
set(gcf,'position',[50,50,1000,800])
plot_meshb_lm_disk(map13,f1_filled,lm1,ShapeIndexHK1);
%title([lr{i},' simulated ferret brain'])
view([90 90])  
caxis([-1,1])
print([name,'_mapping_simu_lm_ShapeIndex.tif'],'-dtiff');


figure(6)
set(gcf,'position',[50,50,1000,800])
plot_meshp_lm_disk(map23,f2_filled,lm2,threshold(mean_curv2,0.35))
%title([lr{i},' simulated ferret brain'])
view([90 90])  
print([name,'_mapping_gel_lm.tif'],'-dtiff');
figure(61)
set(gcf,'position',[50,50,1000,800])
plot_meshp_lm_disk(map23,f2_filled,lm2,ShapeIndexHK2);
%title([lr{i},' simulated ferret brain'])
view([90 90])  
caxis([-1,1])
print([name,'_mapping_gel_lm_ShapeIndex.tif'],'-dtiff');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Evaluate the similarity---shape index
% 03/10/2024
% interpolate the curvature
F = TriScatteredInterp(s3(:,1),s3(:,2),ShapeIndexHK3,'nearest');
interp_ShapeIndexHK3 = F(map13(:,1),map13(:,2));

% perform a binary thresholding
% I1 = (mean_curv1 < mean(mean_curv1)-std(mean_curv1));
% I3 = (interp_mean_curv3 < mean(interp_mean_curv3)-std(interp_mean_curv3));
I1 = ShapeIndexHK1;
I3 = interp_ShapeIndexHK3;

% similarity score
% 2/4/2024 modified--p-norm
score13_p1 = 1 - 0.5*sum(abs(I1-I3))/length(I1)
%score13_p2my = 1 - sqrt( sum( abs(I1-I3).^2 ))/length(I1)
score13_p2 = 1 - 0.5*sqrt( sum(abs(I1-I3).^2)/length(I1) )
score13_p3 = 1 - 0.5*( sum(abs(I1-I3).^3)/length(I1)).^(1/3)
score13_p3 = 1 - 0.5*( sum(abs(I1-I3).^4)/length(I1)).^(1/4)
score13_pinf = 1 - 0.5*max(abs(I1-I3)/length(I1))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F = TriScatteredInterp(s3(:,1),s3(:,2),ShapeIndexHK3,'nearest');
interp_ShapeIndexHK3 = F(map23(:,1),map23(:,2));

% perform a binary thresholding
% I2 = (mean_curv2 < mean(mean_curv2)-std(mean_curv2));
% I3 = (interp_mean_curv3 < mean(interp_mean_curv3)-std(interp_mean_curv3));
I2 = ShapeIndexHK2;
I3 = interp_ShapeIndexHK3;


% similarity score
score23_p1 = 1 - 0.5*sum(abs(I2-I3))/length(I2)
score23_p2 = 1 - 0.5*sqrt( sum(abs(I2-I3).^2)/length(I2) )
score23_p3 = 1 - 0.5*( sum(abs(I2-I3).^3)/length(I2) ).^(1/3)
score23_pinf = 1 - 0.5*max(abs(I2-I3)/length(I2))



fid = fopen('Similarity_ShapeIndex_macaque_Three0426.txt','a+')
fprintf(fid,name)
fprintf(fid,'\r\n')
fprintf(fid,'\t')
fprintf(fid,'simu_p1=%8.4f, gel_p1=%8.4f\r\n',score13_p1,score23_p1)
fprintf(fid,'\t')
fprintf(fid,'simu_p2=%8.4f, gel_p2=%8.4f\r\n',score13_p2,score23_p2)
fprintf(fid,'\t')
fprintf(fid,'simu_p3=%8.4f, gel_p3=%8.4f\r\n',score13_p3,score23_p3)
fprintf(fid,'\t')
% fprintf(fid,'simu_p2my=%8.4f, gel_p2my=%8.4f\r\n',score13_p2my,score23_p2my)
% fprintf(fid,'\t')
fprintf(fid,'simu_pinf=%8.4f, gel_pinf=%8.4f\r\n',score13_pinf,score23_pinf)
fclose('all')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Evaluate the similarity--rescaled mean curvature
% 03/10/2024
% interpolate the curvature
F = TriScatteredInterp(s3(:,1),s3(:,2),mean_curv3,'nearest');
interp_mean_curv3 = F(map13(:,1),map13(:,2));

% perform a binary thresholding
I1 = threshold(mean_curv1,1);
I3 = threshold(interp_mean_curv3,1);

% rescale both curvatures to [0,1]
I1 = (I1-min(I1))/(max(I1)-min(I1))*2-1;
I3 = (I3-min(I3))/(max(I3)-min(I3))*2-1;

% get the interpolated curvature based on the mapping result
% note: the extra factor 0.999 is to avoid interpolation issue at the boundary
% I3_reg = griddata(s2(:,1), s2(:,2), I2, map(:,1)*0.999, map(:,2)*0.999);
% I3_reg(isnan(I2_reg)) = 0;
% 
% % compute the curvature difference on the unit disk
% curv_diff = I1-I2_reg;
% curv_diff
% similarity score
% 2/4/2024 modified--p-norm
score13_p1 = 1 - 0.5*sum(abs(I1-I3))/length(I1)
%score13_p2my = 1 - sqrt( sum( abs(I1-I3).^2 ))/length(I1)
score13_p2 = 1 - 0.5*sqrt( sum(abs(I1-I3).^2)/length(I1) )
score13_p3 = 1 - 0.5*(sum(abs(I1-I3).^3)/length(I1)).^(1/3)

score13_pinf = 1 - 0.5*max(abs(I1-I3)/length(I1))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F = TriScatteredInterp(s3(:,1),s3(:,2),mean_curv3,'nearest');
interp_mean_curv3 = F(map23(:,1),map23(:,2));

% perform a binary thresholding
% I2 = (mean_curv2 < mean(mean_curv2)-std(mean_curv2));
% I3 = (interp_mean_curv3 < mean(interp_mean_curv3)-std(interp_mean_curv3));
% perform a binary thresholding
I2 = threshold(mean_curv2,1);
I3 = threshold(interp_mean_curv3,1);

% rescale both curvatures to [0,1]
I2 = (I2-min(I2))/(max(I2)-min(I2))*2-1;
I3 = (I3-min(I3))/(max(I3)-min(I3))*2-1;


% similarity score
score23_p1 = 1 - 0.5*sum(abs(I2-I3))/length(I2)
score23_p2 = 1 - 0.5*sqrt( sum(abs(I2-I3).^2)/length(I2) )
score23_p3 = 1 - 0.5*(sum(abs(I2-I3).^3)/length(I2) ).^(1/3)

score23_pinf = 1 - 0.5*max(abs(I2-I3)/length(I2))



fid = fopen('Similarity_rescaled_meancurv_macaque_Three0426.txt','a+')
fprintf(fid,name)
fprintf(fid,'\r\n')
fprintf(fid,'\t')
fprintf(fid,'simu_p1=%8.4f, gel_p1=%8.4f\r\n',score13_p1,score23_p1)
fprintf(fid,'\t')
fprintf(fid,'simu_p2=%8.4f, gel_p2=%8.4f\r\n',score13_p2,score23_p2)
fprintf(fid,'\t')
fprintf(fid,'simu_p3=%8.4f, gel_p3=%8.4f\r\n',score13_p3,score23_p3)
fprintf(fid,'\t')
% fprintf(fid,'simu_p2my=%8.4f, gel_p2my=%8.4f\r\n',score13_p2my,score23_p2my)
% fprintf(fid,'\t')
fprintf(fid,'simu_pinf=%8.4f, gel_pinf=%8.4f\r\n',score13_pinf,score23_pinf)
fclose('all')















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Evaluate the similarity
% 03/04/2024
% interpolate the curvature
F = TriScatteredInterp(s3(:,1),s3(:,2),mean_curv3,'nearest');
interp_mean_curv3 = F(map13(:,1),map13(:,2));

% perform a binary thresholding
I1 = (mean_curv1 < mean(mean_curv1)-std(mean_curv1));
I3 = (interp_mean_curv3 < mean(interp_mean_curv3)-std(interp_mean_curv3));

% similarity score
% 2/4/2024 modified--p-norm
score13_p1 = 1 - sum(abs(I1-I3))/length(I1)
score13_p2my = 1 - sqrt( sum( abs(I1-I3).^2 ))/length(I1)
score13_p2 = 1 - sqrt( sum(abs(I1-I3).^2)/length(I1) )
score13_pinf = 1 - max(abs(I1-I3)/length(I1))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 

%% 4/26/2024 adding...
%  Comparison between each two surfaces:
[map12,s1,s2,f1_filled,f2_filled] = disk_landmark_matching_qc_map(v1,f1,lm1,v2,f2,lm2);

% 
figure(90)
set(gcf,'position',[50,50,1000,800])
plot_meshp_lm_disk(s2,f2_filled,lm2,ShapeIndexHK2);
title([lr{i},' simulated human brain12'])
view([90 90])  
%print([name,'_mapping_gel_lm.tif'],'-dtiff');

figure(91)
set(gcf,'position',[50,50,1000,800])
plot_meshg_lm_disk(s1,f1_filled,lm1,ShapeIndexHK1);
title([lr{i},' gel human brain12'])
view([90 90])  
%print([name,'_mapping_real_lm_ShapeIndex.tif'],'-dtiff');


%% Evaluate the similarity---shape index
% 03/10/2024
% interpolate the curvature
F = TriScatteredInterp(s2(:,1),s2(:,2),ShapeIndexHK2,'nearest');
interp_ShapeIndexHK2 = F(map12(:,1),map12(:,2));

I1 = ShapeIndexHK1;
I2 = interp_ShapeIndexHK2;


% similarity score
score12_p1= 1 - 0.5*sum(abs(I1-I2))/length(I1)
score12_p2 = 1 - 0.5*sqrt( sum(abs(I1-I2).^2)/length(I1) )
score12_p3 = 1 - 0.5*( sum(abs(I1-I2).^3)/length(I1)).^(1/3)




fid = fopen('Similarity_ShapeIndexTwo_0426.txt','a+')
fprintf(fid,name)
fprintf(fid,'\r\n')
fprintf(fid,'\t')
fprintf(fid,'simu_gel_p1=%8.4f\r\n',score12_p1)
fprintf(fid,'\t')
fprintf(fid,'simu_gel_p2=%8.4f\r\n',score12_p2)
fprintf(fid,'\t')
fprintf(fid,'simu_gel_p3=%8.4f\r\n',score12_p3)
fprintf(fid,'\t')
fclose('all')





%% 4/26/2024 adding...
%  Comparison between each two surfaces:
[map21,s2,s1,f2_filled,f1_filled] = disk_landmark_matching_qc_map(v2,f2,lm2,v1,f1,lm1);

% 
figure(90)
set(gcf,'position',[50,50,1000,800])
plot_meshp_lm_disk(s2,f2_filled,lm2,ShapeIndexHK2);
title([lr{i},' simulated macaque brain12'])
view([90 90])  
%print([name,'_mapping_gel_lm.tif'],'-dtiff');

figure(91)
set(gcf,'position',[50,50,1000,800])
plot_meshg_lm_disk(s1,f1_filled,lm1,ShapeIndexHK1);
title([lr{i},' gel macaque brain12'])
view([90 90])  
%print([name,'_mapping_real_lm_ShapeIndex.tif'],'-dtiff');


%% Evaluate the similarity---shape index
% 03/10/2024
% interpolate the curvature
F = TriScatteredInterp(s1(:,1),s1(:,2),ShapeIndexHK1,'nearest');
interp_ShapeIndexHK1 = F(map21(:,1),map21(:,2));

I1 = interp_ShapeIndexHK1;
I2 = ShapeIndexHK2;


% similarity score
score21_p1= 1 - 0.5*sum(abs(I1-I2))/length(I1)
score21_p2 = 1 - 0.5*sqrt( sum(abs(I1-I2).^2)/length(I1))
score21_p3 = 1 - 0.5*( sum(abs(I1-I2).^3)/length(I1)).^(1/3)


fid = fopen('Similarity_ShapeIndexTwo_0426.txt','a+')
fprintf(fid,name)
fprintf(fid,'\r\n')
fprintf(fid,'\t')
fprintf(fid,'gel_simu_p1=%8.4f\r\n',score21_p1)
fprintf(fid,'\t')
fprintf(fid,'gel_simu_p2=%8.4f\r\n',score21_p2)
fprintf(fid,'\t')
fprintf(fid,'gel_simu_p3=%8.4f\r\n',score21_p3)
fprintf(fid,'\t')
fclose('all')







