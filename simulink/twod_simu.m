clear
clc

load('Gamma_matrix3.mat')
load('pret2.mat')
let = zeros(21);
for i = 1:21
    let(i) = Gamma_matrix(11,11);
end
pret_2 = pret(:,1);
let_2 = let(:,1);
figure


plot([0:0.05:1],pret_2,'--r','LineWidth',1.2);
hold on;
plot([0:0.05:1],let_2,'-b','LineWidth',1.2);
axis([0 1 1.71 1.764]) 
latexStr1 = ['$\lambda_{11}$ '];
xlabel(latexStr1,'interpreter','latex','FontSize',13) 
latexStr2 = ['$\lambda_{22}$ '];
latexStr3 = ['$\gamma^{*}$ '];
ylabel(latexStr3,'interpreter','latex','FontSize',13) 

latexStr1 = ['$\lambda_{11}=\lambda_{22}$'];                 % LaTeX语法，因为有些东西不用latex写不出来
latexStr2 = ['$\lambda_{11}+\lambda_{22}=1$'];
lgh=legend(latexStr1,latexStr2);  %加注释
set(lgh,'interpreter','latex','FontSize',14)              %就是把字体改成latex格式figure