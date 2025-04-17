% Compare real and simulated brains using disk conformal mapping parameterization
%           
%                              Sifan Yin 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear all; close all;set(0,'defaultfigurecolor','w')

addpath(genpath('Ferret half'));
parentDirectory = fileparts(pwd);
mpath = [parentDirectory,'/FLASH'];
addpath(mpath)

i = 1;
lr = {'left','right'};


name = ['Ferret3_P16_',lr{i}];new_mat = [name,'.mat'];


load(new_mat);
v3 = vertex_real;
v1 = vertex_simu;
v2 = vertex_gel;

f3= face_real;
f1 = face_simu;
f2 = face_gel;



% curvature1 = tricurv(f1,v1);
% mean_curv1 = curvature1.km;
% Gauss_curv1 = curvature1.kg;
% ShapeIndexHK1 = 2/pi*atan(mean_curv1./sqrt(mean_curv1.^2-Gauss_curv1));
% % k11 = curvature1.k1;
% % k12 = curvature1.k2;
% % ShapeIndex1 = 2/pi*atan((k11+k12)./(k12-k11));
% 
% 
% 
% curvature2 = tricurv(f2,v2);
% mean_curv2 = curvature2.km;
% Gauss_curv2 = curvature2.kg;
% ShapeIndexHK2 = 2/pi*atan(mean_curv2./sqrt(mean_curv2.^2-Gauss_curv2));
% k21 = curvature2.k1;
% k22 = curvature2.k2;
% ShapeIndex2 = 2/pi*atan((k21+k22)./(k22-k21));



curvature3= tricurv(f3,v3);
mean_curv3 = curvature3.km;
Gauss_curv3 = curvature3.kg;
ShapeIndexHK3 = 2/pi*atan(mean_curv3./sqrt(mean_curv3.^2-Gauss_curv3));
% k31 = curvature3.k1;
% k32 = curvature3.k2;
% ShapeIndex3 = 2/pi*atan((k31+k32)./(k32-k31));


figure(12)
set(gcf,'position',[50,50,750,600])
%plot_meshg_lm(v3,f3,lm3,ShapeIndex3);
plot_meshj(v3,f3,ShapeIndexHK3);
%title([lr{i},' real ferret brain'])
view([0 10]); caxis([-1,1]) 
% c = colorbar;set(c,'YTick',[-1:0.5:1])
print([name,'_real_lm_SI_jet_left0701.tif'],'-dtiff');



I3 = threshold(mean_curv3,1);
I3r = 2*(I3-min(I3))/(max(I3)-min(I3))-1;

figure(13)
set(gcf,'position',[50,50,750,600])
%plot_meshg_lm(v3,f3,lm3,ShapeIndex3);
plot_meshg(v3,f3,I3r);
%title([lr{i},' real ferret brain'])
view([0 10]); caxis([-1,1]) 
% c = colorbar;set(c,'YTick',[-1:0.5:1])
print([name,'_real_lm_MeanCurv_gray_left0701.tif'],'-dtiff');













% figure(11) % compare mean_curv and shape index
% 
% h3_mean = histogram(threshold0(mean_curv3,10));
% h3_mean.Normalization = 'probability';
% h3_mean.FaceColor = [0 0 0];
% h3_mean.BinWidth = 0.25;
% xlim([-30,max(mean_curv3)])
% title('Mean curvature distribution')
% print('Mean curvature distribution.tif','-dtiff')
% 
% % figure(12)
% % h3_shape = histogram(ShapeIndex3)
% % h3_shape.FaceColor = [0 0 0];
% % h3_shape.Normalization = 'probability';
% % h3_shape.BinWidth = 0.01;
% % title('Shape Index distribution')
% % print('Shape Index distribution.tif','-dtiff')
% 
% figure(13)
% binWidth = 0.005
% Edges = -1:binWidth:1+binWidth;
% h3_shape = histogram(ShapeIndexHK3,Edges)
% h3_shape.FaceColor = [0 0 0];
% h3_shape.Normalization = 'probability';
% h3_shape.BinWidth = binWidth;hold on
% title('Shape index distribution')
% print('Shape Index distribution.tif','-dtiff')
% 


%%%%%%%%%%%%%%%%%%%%%%
%% fitting dist curves
% num_vertex  = length(ShapeIndexHK3);
% counts = histcounts(ShapeIndexHK3,Edges)/num_vertex;
% centers = Edges(1:end-1)+binWidth/2;
% plot(centers,counts,'k.')

%%%%%%%%%%%%%%%%%%%
% figure(10)
% set(gcf,'position',[50,50,750,600])
% binWidth = 0.005
% Edges = -1:binWidth:1+binWidth;
% num_vertex1  = length(ShapeIndexHK1);
% num_vertex2  = length(ShapeIndexHK2);
% num_vertex3  = length(ShapeIndexHK3);
% centers = Edges(1:end-1)+binWidth/2;
% 
% counts1 = histcounts(ShapeIndexHK1,Edges)/num_vertex1;
% counts2 = histcounts(ShapeIndexHK2,Edges)/num_vertex2;
% counts3 = histcounts(ShapeIndexHK3,Edges)/num_vertex3;
% 
% plot(centers,counts3,'k.',centers,counts1,'b.',centers,counts2,'r.','MarkerSize',15);
% xlabel('Shape Index')
% ylabel('Probability')
% xlim([-1 1])
% ylim([0,0.01])
% yticks([0,0.005,0.01])
% %title('Ferret')
% %legend('real','simulated','gel','Location','northwest');
% 
% set(gca,'Fontsize',24,'Fontname','Arial')
% 
% set(0,'defaultlineLinewidth',3);
% set(0,'defaultaxesLinewidth',3);
% set(0,'defaultaxesfontsize',22);
% print('Compare Three shape index_ferret_left_scatter0701.tif','-dtiff')
% 
% %%
% figure(11)
% set(gcf,'position',[50,50,750,600])
% I1 = threshold(mean_curv1,1);
% I2 = threshold(mean_curv2,1);
% I3 = threshold(mean_curv3,1);
% % rescale both curvatures to [0,1]
% I1r = 2*(I1-min(I1))/(max(I1)-min(I1))-1;
% I2r = 2*(I2-min(I2))/(max(I2)-min(I2))-1;
% I3r = 2*(I3-min(I3))/(max(I3)-min(I3))-1;
% 
% binWidth = 0.005
% Edges = -1:binWidth:1+binWidth;
% num_vertex1  = length(ShapeIndexHK1);
% num_vertex2  = length(ShapeIndexHK2);
% num_vertex3  = length(ShapeIndexHK3);
% centers = Edges(1:end-1)+binWidth/2;
% 
% counts1 = histcounts(I1r,Edges)/num_vertex1;
% counts2 = histcounts(I2r,Edges)/num_vertex2;
% counts3 = histcounts(I3r,Edges)/num_vertex3;
% 
% plot(centers,counts3,'k.',centers,counts1,'b.',centers,counts2,'r.','MarkerSize',15);
% xlim([-1,1])
% ylim([0,0.01])
% yticks([0,0.005,0.01])
% xlabel('Rescaled mean curvature')
% ylabel('Probability')
% %title('Mean curvature distribution')
% %legend('real','simulated','gel','Location','northwest');
% set(gca,'Fontsize',24,'Fontname','Arial')
% set(0,'defaultlineLinewidth',3);
% set(0,'defaultaxesLinewidth',3);
% set(0,'defaultaxesfontsize',22);
% print('Compare Three meancurv_ferret_left_scatter0701.tif','-dtiff')
% 
% 
% 
% 
% 
% 
% 
% %% Real ferret jet plot
% figure(12)
% set(gcf,'position',[50,50,750,600])
% %plot_meshg_lm(v3,f3,lm3,ShapeIndex3);
% plot_meshj(v3,f3,ShapeIndexHK3);
% %title([lr{i},' real ferret brain'])
% view([0 10])  
% %colorbar
% print([name,'_real_lm_SI_jet_left0701.tif'],'-dtiff');
% 
% 

function v_chopped = threshold(v, cutoff)
v_chopped = v;
v_chopped(v > mean(v) + cutoff*std(v)) = mean(v) + cutoff*std(v);
v_chopped(v < mean(v) - cutoff*std(v)) = mean(v) - cutoff*std(v);
end











