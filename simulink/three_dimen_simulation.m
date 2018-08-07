clear
clc

load("Gamma_matrix.mat")

figure
mesh([0:0.05:1],[0:0.05:1],Gamma_matrix);% toc
%  axis([0 1 0.2 0.3]);
latexStr1 = ['$\lambda_{11}$ '];
xlabel(latexStr1,'interpreter','latex','FontSize',14) 
latexStr2 = ['$\lambda_{22}$ '];
ylabel(latexStr2,'interpreter','latex','FontSize',14) 
latexStr3 = ['$\gamma^*$ '];
zlabel(latexStr3,'interpreter','latex','FontSize',14) 
view(60,10)