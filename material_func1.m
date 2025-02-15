function [mat_piezok1s,mat_flexok1s,mat_piezo1,mat_flexo1,mat_high]=material_func1()
%% Material (2) BaTio3
%----------edit-------------
% on 24-feb-2024
% added post processing h term check paper 2 draft 2 page 6
% %BaTio3 
%elastic
c11 = 16.6e10;     
c13 = 7.75e10;         
c33 = 16.2e10;     
c55 = 4.29e10;

%piezo
e31 = -4.4 ;    
e33 = 18.6;    
e15 = 11.6; 

% dielectric
di11 = 1.2558e-08 ;%8.85e-12*1419 ;   % Dielectric constants 8.85e-12*1268, 1500-2000 range
di33 = 1.1231e-08 ;%8.85e-12*1269 ;

%flexoelectric:

nu11=0;   %original 50e-6 
nu12=1e-6;    % Paper3=3e-6 1e-6%1.4e-7  %1.4e-6,2e-6   small value 1.4e-8 is for field plot original 50e-6  use 10e-6
nu44=0;

%%

elastic = [c11 c13 0; c13 c33 0; 0 0 c55];

% length scale

% for PZT5h internal length parameter -Sladek-2017
%qgr=1;
%lh=sqrt(qgr)*5*10e-9;
%Sizf=lh^2;


% legth scale Ghasemi-2019

lh=5e-9;     
Sizf=lh^2;    
h_mat1=zeros(6,6);
h_mat1(1,1)=c11;h_mat1(2,2)=c11;h_mat1(3,3)=c55;h_mat1(4,4)=c11;h_mat1(5,5)=c11;h_mat1(6,6)=c55;
h_mat1(1,2)=c13;h_mat1(2,1)=c13;h_mat1(4,5)=c13;h_mat1(5,4)=c13;

h_mat=Sizf*h_mat1; %  

%BaTio3

piezo = [0 e31
         0 e33
         e15 0] ;
% conditioning scaling Refer paper2_Draft2
betcon=10e10;
piezos = betcon*piezo;   
%BaTio3
  
di = [(di11) 0 ; 0 (di33)] ;

% conditioning scaling Refer paper2_Draft2

dis = (betcon^2)*di;   


%flexo 

nu_s=[nu11 nu12 0 0 0 nu44
      0 0 nu44 nu12 nu11 0];

nu_ss = betcon*nu_s; 



 
%material matrices


% materials without conditioning

% mat_piezok2= [elastic piezo;
%                piezo' -di];        % for K matrix
% mat_flexok2=[h_mat nu_s';          % replace zeros(6,6) by h_mat if u want to add scale vversa
%              nu_s zeros(2,2)];

         
 % materials with conditioning
mat_piezok1s= [elastic piezos
               piezos' -dis];        % for K matrix
mat_flexok1s=[h_mat nu_ss'          % replace zeros(6,6) by h_mat if u want to add scale vversa
             nu_ss zeros(2,2)];
 


             
         
%% post processing material constants material (2)
% no scaling is required here

he=zeros(5,12);
he(1,1)=c11;he(1,2)=c13;he(1,10)=c11;he(1,11)=c13;
he(2,1)=c13;he(2,2)=c33;he(2,10)=c13;he(2,11)=c33;
he(3,3)=c55;he(3,12)=c55;
hel=Sizf*he;

nu_e=[nu11 0 0 nu12
    nu12 0 0 nu11
    0 nu44 nu44 0];


mat_piezo1 = [elastic piezo
               piezo' -di] ;        %for post processing(this is verified)
mat_flexo1=[zeros(3,6) -nu_e
             nu_s zeros(2,4)];
%% for Higher order stress sigma_ijk

mat_high=[h_mat nu_s'];