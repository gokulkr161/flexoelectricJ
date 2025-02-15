%---------------------------------------------
% Compute J integral Flexo problem of higher order
% this function is for center FGM problem 
function [JintMech,JintPiez,JintFlex]=JintCentFlexoFunc_FGMhighO(Jdomain,qnode,QT,elementV,split_elem,uKnot,vKnot,enrich_node,controlPts,element,...
    p,q,weights,elRangeU,elRangeV,node,index,pos,U,D,xtip,xCr)

% I1 = 0;
% I2 = 0;
% I  = [zeros(2,1)];
%Jint=0;ind=0;ind1=0;

Ideb=[];
JintFlex=0;
JintMech=0;
JintPiez=0;
globGP=[];

% ---------------------------
% Starting LOOP over ELEMENTS
%----------------------------

for iel = 1 : size(Jdomain,2)           % no of elements in J domain
    e      = Jdomain(iel) ;    % current element
    sctr   = element(e,:);     % 9 control points of IGA element 
    sctrQ4 = elementV(e,:);    % 4 nodes of postprocessing mesh
    nn     = length(sctr);      % 9
    
    idu    = index(e,1);
    idv    = index(e,2);
    xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]   Knot span or IGA element span
    etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
    
    % Choose Gauss quadrature rule
    
    if (ismember(e,split_elem))     % split element
        order = 12; % 13 GPs for each sub-triangle
        %phi   = ls(sctr,1);
        %[W,Q] = discontQ4quad(order,phi);              % For triangulation  Delaunay
        [W,Q] = quadrature(order,'GAUSS',2);
    else
        order = 8 ; 
        [W,Q] = quadrature(order,'GAUSS',2);
    end
    
    % -----------------------------
    % start loop over Gauss points
    % -----------------------------
    
    for igp = 1:size(W,1)        
        pt = Q(igp,:);
        wt = W(igp);
        
        % Q4 element for weighting q   ** We are using lagrange for
        % interpolate weighting function q
        
        [N,dNdxi] = lagrange_basis('Q4',pt);
        J0    = node(sctrQ4,:)'*dNdxi;     % Jacobian for 4 noded quadrilateral element
        invJ0 = inv(J0);
        dNdx  = dNdxi*invJ0;
        
        % compute coords in parameter space
        Xi      = parent2ParametricSpace(xiE,pt(1));
        Eta     = parent2ParametricSpace(etaE,pt(2));
        J2      = jacobianPaPaMapping(xiE,etaE);
    
        % compute derivative of basis functions w.r.t parameter coord


      [N, dRdxi, dRdeta, dR2dxi, dR2deta, dR2dxideta] = NURBS2DBasis2ndDers( [Xi; Eta], p, q, uKnot,vKnot, weights);
        pts        = controlPts(sctr,:);
        jacob      = pts'*[dRdxi' dRdeta'];
%         J1         = det(jacob);
        invJacob   = inv(jacob);
        dRdx       = [dRdxi' dRdeta'] * invJacob;
        Gpt        = N*pts;     % GP in global coord
        
        globGP     = [globGP; Gpt]; % for plotting GPs only
        
        % +++++++++++++++++++++++++ 
        % Gradient of displacement
        % +++++++++++++++++++++++++ 
        
        % need to compute u,x u,y v,x v,y, stored in matrix H
        
        [B_piezo,B_flexo,B_flexoPost,J1] = BForJintCenter(Xi,Eta,e,enrich_node,xCr,xtip,N,dRdxi,dRdeta,dR2dxi,dR2deta,dR2dxideta);
        leB = size(B_piezo,2);
        
        % nodal displacement of current element
        % taken from the total nodal parameters U
        
        elemDisp = element_disp(e,pos,enrich_node,U);
        
        % compute derivatives of u w.r.t x and y
        % gradient of scalar field potential is a vecor field
        
        H(1,1) = B_piezo(1,1:3:leB)*elemDisp(1:3:leB);    % u,x
        H(1,2) = B_piezo(2,2:3:leB)*elemDisp(1:3:leB);    % u,y
        
        H(2,1) = B_piezo(1,1:3:leB)*elemDisp(2:3:leB);    % v,x
        H(2,2) = B_piezo(2,2:3:leB)*elemDisp(2:3:leB);    % v,y
        
        dp_dx = B_piezo(1,1:3:leB)*elemDisp(3:3:leB);    % electric potential,x  
        dp_dy = B_piezo(2,2:3:leB)*elemDisp(3:3:leB);    % electric potential,y
        Grad_potential=[dp_dx;dp_dy]; % we need dp_dx only

        % +++++++++++++++++++
        % Second gradient of displacement
        % +++++++++++++++++++
        HH(1,1) = B_flexo(1,1:3:leB)*elemDisp(1:3:leB);    % u,xx
        HH(1,2) = B_flexo(2,2:3:leB)*elemDisp(1:3:leB);    % u,yx

        HH(2,1) = B_flexo(1,1:3:leB)*elemDisp(2:3:leB);    % v,xx
        HH(2,2) = B_flexo(2,2:3:leB)*elemDisp(2:3:leB);    % v,yx
        
        % making the 3rd order tesnor complete
        HH(:,:,2)=zeros(2,2);   % initiation of second array of 3rd order tensor

        HH(1,1,2) = B_flexo(4,1:3:leB)*elemDisp(1:3:leB);    % u,xy
        HH(1,2,2) = B_flexo(5,2:3:leB)*elemDisp(1:3:leB);    % u,yy

        HH(2,1,2) = B_flexo(4,1:3:leB)*elemDisp(2:3:leB);    % v,xy
        HH(2,2,2) = B_flexo(5,2:3:leB)*elemDisp(2:3:leB);    % v,yy



        % +++++++++++++++++++
        % Gradient of weight
        % +++++++++++++++++++ 
        
        weight  = qnode(iel,:);   % q is the weighting function 
        %gradq   = weight*dRdx;
        gradq   = weight*dNdx;

        
        
%         % +++++++++++++++++++
%         %        FGM
%         % +++++++++++++++++++ 
        % material_3;   PZT 5h
        [~,~,mat_piezo1,mat_flexo1,mat_high]=material_func1(); %BT
        % material_grad_PZT  % only for validation
%           Gausspoint = N *pts;
%           
%        YP=Gausspoint(2);
%        [mat_piezok1s,mat_flexok1s,mat_piezo1,mat_flexo1]=material_func1()
%       [mat_piezo1,mat_flexo1]=FGM(mat_piezo1,mat_flexo1,mat_piezo2,mat_flexo2,D,YP);
% 
% 
%         
%         
        
        
        
        

        % ++++++++++++++
        % Stress at GPs
        % ++++++++++++++ 
                    strainPiezo          = B_piezo*elemDisp;
            strainGrad_Flexo     = B_flexoPost*elemDisp;    % H matrix* Disp      % stress and strain at nodes
            sigma           = (mat_piezo1*strainPiezo + mat_flexo1*strainGrad_Flexo);  % flexo piezo coupled
            electrDisp  = sigma(4:5);  % electric displacement
            sigma= sigma(1:3);   % mechanical stresses
            epsilon=strainPiezo(1:3);
            elecField=strainPiezo(4:5);
        %         stress(e,gp,1:5)= sigma;      % 3 stresses and 2 electric displacements
        %         epsilon = B*elemDisp ;
        %         sigma   = C*epsilon;

        % higher order stress
        Sgrad_E=B_flexo*elemDisp;       % vector containing strain grad and E
        sigmahi=mat_high*Sgrad_E;       % third order stress sigma_ijk
        epsilon_grad= strainGrad_Flexo(1:6);  % strain grad for strain energy dens calculation
        %strain_grad=Sgrad_E(1:6);   % this also gives  strain grad

        
        % +++++++++++++++++++++++++++++++++++
        % Transformation to local coordinate Stress transformation
        % +++++++++++++++++++++++++++++++++++ 
        
        voit2ind    = [1 3;3 2];                  % to convert stress vector to stress matrix
        gradqloc    = QT*gradq';                             
        graddisploc = QT*H*QT';                   % gradient of disp field
        GradPotloc=QT*Grad_potential;            % gradient of potential                                    
        stressloc   = QT*sigma(voit2ind)*QT';      
        elecDisploc= QT*electrDisp;
        epsilon(3)  = 0.5 * epsilon(3);    % in strain energy dencity we are using .5 engineerng strain
        strainloc   = QT*epsilon(voit2ind)*QT';
        elecFieldploc=QT*elecField;

        % +++++++++++++++++++++++++++++++++++
        % Transformation to local coordinate of higher order tensors
        % +++++++++++++++++++++++++++++++++++
        voit2ind3=voit2ind;   %initiation
        voit2ind3(:,:,2)=[4 6;6 5];   % position 3rd order tensor ie two 2x2 matrieces in matlab
        sigmahitens=sigmahi(voit2ind3);    % higher order stress in 2x2x2 form
        epsilon_grad(3)=0.5*epsilon_grad(3);
        epsilon_grad(6)=0.5*epsilon_grad(6);
        epsilon_gradtens=epsilon_grad(voit2ind3); %  strain grads in 2x2x2 form
        sigmahiloc=thirdTranf(sigmahitens,QT);     % sigmahi transformed to local
        epsilon_gradloc=thirdTranf(epsilon_gradtens,QT);  % strain grad transformed to local
        gragradisp_loc=thirdTranf(HH,QT);               % second gradient of disp is tranformed

        % +++++++++++++++
        %   J integral
        % +++++++++++++++
        % mechanical energy
        dSEnmech = 0;
        for i=1:2
            for j=1:2
                dSEnmech = dSEnmech + 0.5*(stressloc(i,j)*strainloc(i,j)); 
            end
        end
        % electrical energy density
        dSEnelec=0;
        for i=1:2
            dSEnelec = dSEnelec+0.5*(elecDisploc(i)*elecFieldploc(i)); 
        end
        % higher order energy density
        dSEnhigh=0;
        for i=1:2
            for j=1:2
                for k=1:2
                    dSEnhigh = dSEnhigh + 0.5*(sigmahiloc(i,j,k)*epsilon_gradloc(i,j,k)); 
                end
            end
        end

        dSEnpiezo=dSEnmech-dSEnelec;         % pure piezo strain energy
        dSEnergy=dSEnmech+dSEnhigh-dSEnelec; % total strain energy density
        %Mech part of J-integral
        JJmech= (stressloc(1,1) * graddisploc(1,1) + stressloc(2,1) * graddisploc(2,1)) * gradqloc(1) + ...
            (stressloc(1,2) * graddisploc(1,1) + stressloc(2,2) * graddisploc(2,1)) * gradqloc(2);
        %piezo part of J-integral

        JJp= ((stressloc(1,1) * graddisploc(1,1) + stressloc(2,1) * graddisploc(2,1))+elecDisploc(1)*GradPotloc(1) ) * gradqloc(1) + ...
            ((stressloc(1,2) * graddisploc(1,1) + stressloc(2,2) * graddisploc(2,1))+elecDisploc(2)*GradPotloc(1) ) * gradqloc(2);
        
        %higher order part of J-integral

        JJh=0;
        for i=1:2
            for j=1:2
                for k=1:2
                    JJh = JJh + sigmahiloc(i,k,j)*gragradisp_loc(i,k,1)*gradqloc(j); 
                end
            end
        end
        JJ1=JJp+JJh;               % J flexo
        % J piezo J flexo J-mech
        JintFlex = JintFlex + ((JJ1-dSEnergy*gradqloc(1))*J1*J2*wt);      % integration for flexo J
        JintPiez = JintPiez + ((JJp-dSEnpiezo*gradqloc(1))*J1*J2*wt);      % integration for piezo J
        JintMech = JintMech + ((JJmech-dSEnmech*gradqloc(1))*J1*J2*wt);      % integration for mech J
        
         %--------------%--------------------%--------------------%---------------------%       
          

    end       % of quadrature loop
end           % end of element loop


%% plot
% figure
% hold on
% % plot the circle
% theta = -pi:0.1:pi;
% xo = xtip(1) + radius*cos(theta) ;
% yo = xtip(2) + radius*sin(theta) ;
% plot(xo,yo,'k-');
% plot_mesh(node,elementV,'Q4','b-',1.2)
% plot_mesh(node,elementV(Jdomain,:),'Q4','r-',1.2)
% cr = plot(xCr(:,1),xCr(:,2),'k-');
% set(cr,'LineWidth',2);
% cr = plot(globGP(:,1),globGP(:,2),'b*');
% -------------------------------------