% Ai_1 = [-0.4799 5.1546; -3.8162 14.6732]
% Ai_2 = [-1.6026 9.1632; -0.5918 3.0317]
% Ai_3 = [0.6346 0.9178; -0.5056 2.4811]

clear;
clc;

size_x = 2;
size_y = 1;
size_u = 1;
size_z = 1;

N = 2;  % number of mode i.
M = 2;  % number of mode m.

A1 = [0.4 0.40; 0.2 1];% A2 = [1.58 -0.36; 0.80 -1.12];
A2 = [1.1 0.6; 0.3 0.4];
A = cell(1,N); A{1,1} = A1;  A{1,2} = A2;

B1 = [0.5; 0.5]; B2 = [0.7; 0.5];
B = cell(1,N); B{1,1} = B1; B{1,2} = B2;

F1 = [1; 1.2]; F2 = [1.2; 1]; 
F = cell(1,N); F{1,1} = F1; F{1,2} = F2;

C1 = [0.9 0.5]; C2 = [-1 0.7]; 
C = cell(1,N); C{1,1} = C1; C{1,2} = C2;

C1z = [1 1]; C2z = [1 1];
Cz = cell(1,N); Cz{1,1} = C1z; Cz{1,2} = C2z;

E1x = [0.1; 0.12]; E2x = [0.1; 0.11]; 
Ex = cell(1,N); Ex{1,1} = E1x; Ex{1,2} = E2x;

iOmeg1 = 1.3;
iOmeg2 = 1.5;
IOmega = cell(1,N);
IOmega{1,1} = iOmeg1;
IOmega{1,2} = iOmeg2;

P_Matrix = [0.6 0.4; 0.2 0.8];
% M_Matrix = [0.6 0.4; 0.7 0.3];
% M_Matrix = [0.4 0.6; 0.5 0.5];
% M_Matrix = [0 1; 1 0];
M_Matrix = [1 0; 0 1];

P = cell(1,N);
R11 = cell(N,M);
R12 = cell(N,M);
R22 = cell(N,M);
T = cell(1,N);
K = cell(1,M);
G = cell(1,M);
K_sp = cell(1,N);
G_sp = cell(1,N);
G_rv = cell(1,M);
K_rv = cell(1,M);
% 开始LMI的部分
setlmis([]);    %以setlmis开始

for i = 1:N    %定义barP_i
    P{1,i} = lmivar(1,[size_x 1]);
    T{1,i} = lmivar(1,[size_y 1]);
end

for m = 1:M
    K{1,m} = lmivar(2,[size_u size_x]);
    G{1,m} = lmivar(2,[size_u size_y]);
end

for i = 1:N
    for m = 1:M
        R11{i,m} = lmivar(1,[size_x 1]);
        R12{i,m} = lmivar(2,[size_x size_y]);
        R22{i,m} = lmivar(1,[size_y 1]);
    end
end


k= 1;
for i = 1:N
    lmiterm([k 1 1 P{1,i}],-1,1);    %pi>0
    k = k+1; 
    lmiterm([k 1 1 T{1,i}],-1,1);    %Ti>0
    k = k+1;
end
for i = 1:N
    for m = 1:M   %Rim>0
        lmiterm([k 1 1 R11{i,m}],-1,1);
        lmiterm([k 1 2 R12{i,m}],-1,1);
        lmiterm([k 2 2 R22{i,m}],-1,1);
        k = k+1;
    end
end

for i = 1:N
    for m = 1:M
        % define condition 1 in the follows.
        lmiterm([k 1 1 R11{i,m}],-1,1);    
        lmiterm([k 1 2 R12{i,m}],-1,1);
        lmiterm([k 2 2 R22{i,m}],-1,1);
        for ii = 1:N
            lmiterm([k 1 ii+2 0],sqrt(P_Matrix(i,ii))*A{1,i}');
            lmiterm([k 2 ii+2 0],sqrt(P_Matrix(i,ii))*F{1,i}');
            lmiterm([k 1 ii+2 -K{1,m}],sqrt(P_Matrix(i,ii)),B{1,i}');
            lmiterm([k 2 ii+2 -G{1,m}],sqrt(P_Matrix(i,ii)),B{1,i}');
            lmiterm([k ii+2 ii+2 P{1,ii}],-1,1);
        end
        k = k+1;
        
        
    end
    
    %define condition 2 as follows.
    lmiterm([k 1 1 P{1,i}],-1,1);
    lmiterm([k 2 2 0],-eye(size_x));
    lmiterm([k 2 3 T{1,i}],C{1,i}'*IOmega{1,i},1);
    lmiterm([k 3 3 T{1,i}],-2,1);
    for mm = 1:M
        lmiterm([k 1 3*mm+1 P{1,i}],sqrt(M_Matrix(i,mm)),1);
        lmiterm([k 2 3*mm+2 R11{i,mm}],sqrt(M_Matrix(i,mm)),1);
        lmiterm([k 2 3*mm+3 R12{i,mm}],sqrt(M_Matrix(i,mm)),1);
        lmiterm([k 3 3*mm+3 R22{i,mm}],sqrt(M_Matrix(i,mm)),1);
        lmiterm([k 3*mm+1 3*mm+1 0],-eye(size_x));
        lmiterm([k 3*mm+2 3*mm+2 R11{i,mm}],-1,1);
        lmiterm([k 3*mm+2 3*mm+3 R12{i,mm}],-1,1);
        lmiterm([k 3*mm+3 3*mm+3 R22{i,mm}],-1,1);
    end
    k = k+1;
end


lmisys = getlmis;
[tmin,xfeas] = feasp(lmisys);
% Ps = dec2mat(lmisys,xfeas,P{1,1})
% ks = dec2mat(lmisys,xfeas,R11{1,1})


        





























