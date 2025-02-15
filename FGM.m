% 


function [mat_piezo1,mat_flexo1]=  FGM(mat_piezo1,mat_flexo1,mat_piezo2,mat_flexo2,h,Y)


n=7;    % change
VFAC=Y/h;
[mat_piezo1]=FGM_avg(mat_piezo1,mat_piezo2,n,VFAC);
[mat_flexo1]=FGM_avg(mat_flexo1,mat_flexo2,n,VFAC);

 function [E]=FGM_avg(E1,E2,n,VFAC)
    E=(E2-E1)*((VFAC)^n)+E1; % here VFAC is the vol fractn
 end
end


% for ENUM=1:T_ENUM
% 
% if n<100
%             COR= N_SHAP * ELE_CORRY';
%             YCOR=COR;
%             ETAB=(2*YCOR-h)/(h);
%             E=(E1-E2)*((1+ETAB)/2)^n+E2;
% else
%             E=E2;
% end
%         C = elasticityMatrix(E,nu,stressState);
% end
