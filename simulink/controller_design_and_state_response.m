clear;
clc;

size_x = 2;
size_y = 1;
size_u = 1;
size_z = 1;
size_w = 1;



A1 = [0.4 0.40; 0.2 1];
A2 = [1.1 0.6; 0.3 0.4];
A{1} = A1;  A{2} = A2;

B1 = [0.5; 0.5]; B2 = [0.7; 0.5];
B{1} = B1; B{2} = B2;

F1 = [1; 1.2]; F2 = [1.2; 1]; 
F{1} = F1; F{2} = F2;

C1 = [0.9 0.5]; C2 = [-1 0.7]; 
C{1} = C1; C{2} = C2;



% C1z = [1 1]; C2z = [1 1];
C1z = [0.2  0]; C2z = [0.1 0.1];
% C1z = [0.3  0.2]; C2z = [0.1 0.2];
% C1z = [0 0]; C2z = [0 0];
Cz{1} = C1z; Cz{2} = C2z;

% E1x = [0.1; 0.12]; E2x = [0.1; 0.11]; 
E1x = [1; 0.5]; E2x = [1; 0.5]; 
% E1x = [0; 0]; E2x = [0; 0]; 
Ex{1} = E1x; Ex{2} = E2x;

Gz{1} = 0.3; Gz{2} = 0.4;
Dz{1} = 0; Dz{2} = 0.3;
Ez{1} = 1.3; Ez{2} = -0.8;

% Gz{1} = 1; Gz{2} = 2;
% Dz{1} = 1; Dz{2} = 2;
% Ez{1} = 0.5; Ez{2} = 0.5;

iOmeg1 = 1.3;
iOmeg2 = 1.5;
IOmega{1} = iOmeg1;
IOmega{2} = iOmeg2;

% N = 2;  % number of mode i.
% M = 2;  % number of mode m.
N = 2;  % number of mode i.
M = 2;  % number of mode m.

P_Matrix = [0.6 0.4; 0.2 0.8];
% M_Matrix = [ 0.4 0.6; 0.3 0.7];
% M_Matrix = [0.4 0.6; 0 1];
M_Matrix = [1 0; 0 1];
% M_Matrix = [1 0 0; 0 1 0; 0 0 1];
% M_Matrix = [0 0 1; 0 1 0; 0 0 1];
% M_Matrix = [0.8 0.1 0.1; 0.2 0.3 0.5; 0 0 1];
% M_Matrix = [0.8 0.1 0.1; 0.2 0.3 0.5; 0.1 0.8 0.1];
% 开始LMI的部分
setlmis([]);    %以setlmis开始

for i = 1:N    %定义barP_i
    P{i} = lmivar(1,[size_x 1]);
    T{i} = lmivar(1,[size_y 1]);
end

for m = 1:M
    K{m} = lmivar(2,[size_u size_x]);
    G{m} = lmivar(2,[size_u size_y]);
end

for i = 1:N
    for m = 1:M
        R11{i,m} = lmivar(1,[size_x 1]);
        R12{i,m} = lmivar(2,[size_x size_y]);
        R22{i,m} = lmivar(1,[size_y 1]);
        R13{i,m} = lmivar(2,[size_x size_w]);
        R23{i,m} = lmivar(2,[size_y size_w]);
        R33{i,m} = lmivar(1,[size_w 1]);
    end
end

varGamma = lmivar(1,[1,0]);


k= 1;
for i = 1:N
    lmiterm([k 1 1 P{i}],-1,1);    %pi>0
    k = k+1; 
    lmiterm([k 1 1 T{i}],-1,1);    %Ti>0
    k = k+1;
end
for i = 1:N
    for m = 1:M   %Rim>0
        lmiterm([k 1 1 R11{i,m}],-1,1);
        lmiterm([k 1 2 R12{i,m}],-1,1);
        lmiterm([k 2 2 R22{i,m}],-1,1);
        lmiterm([k 1 3 R13{i,m}],-1,1);
        lmiterm([k 2 3 R23{i,m}],-1,1);
        lmiterm([k 3 3 R33{i,m}],-1,1);
        k = k+1;
    end
end

lmiterm([k 1 1 varGamma],-1,1);
k = k+1;

for i = 1:N
    for m = 1:M
        % define condition 1 in the follows.
%         lmiterm([k 1 1 R11{i,m}],-1,1);
%         lmiterm([k 1 2 R12{i,m}],-1,1);
%         lmiterm([k 2 2 R22{i,m}],-1,1);
%         for ii = 1:N
%             lmiterm([k 1 ii+2 0],sqrt(P_Matrix(i,ii))*A{i}');
%             lmiterm([k 2 ii+2 0],sqrt(P_Matrix(i,ii))*F{i}');
%             lmiterm([k 1 ii+2 -K{m}],sqrt(P_Matrix(i,ii)),B{i}');
%             lmiterm([k 2 ii+2 -G{m}],sqrt(P_Matrix(i,ii)),B{i}');
%             lmiterm([k ii+2 ii+2 P{ii}],-1,1);
%         end
        
        lmiterm([k 1 1 R11{i,m}],-1,1);
        lmiterm([k 1 2 R12{i,m}],-1,1);
        lmiterm([k 2 2 R22{i,m}],-1,1);
        lmiterm([k 1 3 R13{i,m}],-1,1);
        lmiterm([k 2 3 R23{i,m}],-1,1);
        lmiterm([k 3 3 R33{i,m}],-1,1);
        for ii = 1:N
            lmiterm([k 1 ii+3 0],sqrt(P_Matrix(i,ii))*A{i}');
            lmiterm([k 2 ii+3 0],sqrt(P_Matrix(i,ii))*F{i}');
            lmiterm([k 3 ii+3 0],sqrt(P_Matrix(i,ii))*Ex{i}');
            lmiterm([k 1 ii+3 -K{m}],sqrt(P_Matrix(i,ii)),B{i}');
            lmiterm([k 2 ii+3 -G{m}],sqrt(P_Matrix(i,ii)),B{i}');
            lmiterm([k ii+3 ii+3 P{ii}],-1,1);
        end
        lmiterm([k 1 ii+4 0],Cz{i}');
        lmiterm([k 2 ii+4 0],Gz{i}');
        lmiterm([k 3 ii+4 0],Ez{i}');
        lmiterm([k 1 ii+4 -K{m}],1,Dz{i}');
        lmiterm([k 2 ii+4 -G{m}],1,Dz{i}');
        lmiterm([k ii+4 ii+4 0],-eye(size_w));
        k = k+1;
    end
    
    %define condition 2 as follows.
    lmiterm([k 1 1 P{i}],-1,1);
    lmiterm([k 2 2 0],-eye(size_x));
    lmiterm([k 2 3 T{i}],C{i}'*IOmega{i},1);
    lmiterm([k 3 3 T{i}],-2,1);
    lmiterm([k 4 4 varGamma],-eye(size_w),1);
    for mm = 1:M
%         lmiterm([k 1 3*mm+1 P{i}],sqrt(M_Matrix(i,mm)),1);
%         lmiterm([k 2 3*mm+2 R11{i,mm}],sqrt(M_Matrix(i,mm)),1);
%         lmiterm([k 2 3*mm+3 R12{i,mm}],sqrt(M_Matrix(i,mm)),1);
%         lmiterm([k 3 3*mm+3 R22{i,mm}],sqrt(M_Matrix(i,mm)),1);     
%         lmiterm([k 3*mm+1 3*mm+1 0],-eye(size_x));
%         lmiterm([k 3*mm+2 3*mm+2 R11{i,mm}],-1,1);
%         lmiterm([k 3*mm+2 3*mm+3 R12{i,mm}],-1,1);
%         lmiterm([k 3*mm+3 3*mm+3 R22{i,mm}],-1,1);
        
        lmiterm([k 1 4*mm+1 P{i}],sqrt(M_Matrix(i,mm)),1);
        lmiterm([k 2 4*mm+2 R11{i,mm}],sqrt(M_Matrix(i,mm)),1);
        lmiterm([k 2 4*mm+3 R12{i,mm}],sqrt(M_Matrix(i,mm)),1);
        lmiterm([k 3 4*mm+3 R22{i,mm}],sqrt(M_Matrix(i,mm)),1);
        lmiterm([k 2 4*mm+4 R13{i,mm}],sqrt(M_Matrix(i,mm)),1);
        lmiterm([k 3 4*mm+4 R23{i,mm}],sqrt(M_Matrix(i,mm)),1);
        lmiterm([k 4 4*mm+4 R33{i,mm}],sqrt(M_Matrix(i,mm)),1);
        
        lmiterm([k 4*mm+1 4*mm+1 0],-eye(size_x));
        lmiterm([k 4*mm+2 4*mm+2 R11{i,mm}],-1,1);
        lmiterm([k 4*mm+2 4*mm+3 R12{i,mm}],-1,1);
        lmiterm([k 4*mm+2 4*mm+4 R13{i,mm}],-1,1);
        lmiterm([k 4*mm+3 4*mm+3 R22{i,mm}],-1,1);
        lmiterm([k 4*mm+3 4*mm+4 R23{i,mm}],-1,1);
        lmiterm([k 4*mm+4 4*mm+4 R33{i,mm}],-1,1);
    end
    k = k+1;
end


lmisys = getlmis;
% [tmin,xfeas] = feasp(lmisys);
% Ps = dec2mat(lmisys,xfeas,P{1})
% ks = dec2mat(lmisys,xfeas,R11{1,2})
% gag = dec2mat(lmisys,xfeas,varGamma)

var_count = decnbr(lmisys);
c = [zeros(var_count-1,1);1];
[xopt,bpot] = mincx(lmisys,c,[0 200 0 0 0 ]);

%状态反馈及输出反馈求解，还有性能指标。
result_gamma = dec2mat(lmisys,bpot,varGamma)    
for i=1:M
    controllerK{i} = dec2mat(lmisys,bpot,K{i});
    controllerG{i} = dec2mat(lmisys,bpot,G{i});
end



%%下面要开始做仿真了哟。
%先做一个状态变化吧.
final_time  = 100;
pmode = [];
kmode = [];
pmode(1) = 1;
for i=2:1:final_time
    sk = pmode(i-1);
    if rand<P_Matrix(sk,1)
        pmode(i) = 1;
    else
        pmode(i) = 2;
    end
end
for i=1:1:final_time
    sk = pmode(i);
    if rand<M_Matrix(sk,1)
        kmode(i) = 1;
    else
        kmode(i) = 2;
    end
end

% Generate disturbance input.
w = [];
for i=1:final_time
    w(i) = sin(i)*0.9^i;
end

figure
stairs([1:1:final_time],pmode);
figure
stairs([1:1:final_time],kmode);
% figure
% plot([1:1:final_time],w);

x{1} = [2; -2.5];
u{1}=0;
for i=1:final_time-1
    mode_p = pmode(i);
    mode_k = kmode(i);
    if mode_p==1
        x{i+1} = (A{mode_p}+B{mode_p}*controllerK{mode_k})*x{i}+(F{mode_p}+B{mode_p}*controllerG{mode_k})*0.5*iOmeg1*C{mode_p}*x{i}*(1+cos(25*C{mode_p}*x{i}))+Ex{mode_p}*w(i);
        u{i} = controllerK{mode_k}*x{i}+controllerG{mode_k}*0.5*iOmeg1*C{mode_p}*x{i}*(1+cos(25*C{mode_p}*x{i}));
    else
        x{i+1} = (A{mode_p}+B{mode_p}*controllerK{mode_k})*x{i}+(F{mode_p}+B{mode_p}*controllerG{mode_k})*0.5*iOmeg2*C{mode_p}*x{i}*(1-sin(20*C{mode_p}*x{i}))+Ex{mode_p}*w(i);
        u{i} = controllerK{mode_k}*x{i}+controllerG{mode_k}*0.5*iOmeg2*C{mode_p}*x{i}*(1-sin(20*C{mode_p}*x{i}));
    end
%     x{i+1} = (A{mode_p}+B{mode_p}*controllerK{mode_k})*x{i}+(F{mode_p}+B{mode_p}*controllerG{mode_k})*C{mode_p}*x{i}+Ex{mode_p}*w(i);
end
u{100}=0;
for i=1:final_time
    x1(i) = x{i}(1);
    x2(i) = x{i}(2);
    us(i) = u{i}(1);
end
figure
% plot([1:1:final_time],x1); hold on;
% plot([1:1:final_time],x2); 
stairs([1:1:final_time],x1,'-'); hold on;
stairs([1:1:final_time],x2,'-.');
figure
stairs([1:1:final_time],us,'-');