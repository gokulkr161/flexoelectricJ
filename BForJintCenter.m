function [B_piezo,B_flexo,B_flexoPost,J1] = BForJintCenter(Xi,Eta,e,enrich_node,xCr,xTip,N,dRdxi,dRdeta,dR2dxi,dR2deta,dR2dxideta)
%modified by Gokul IIT Mandi 10-2021

% this code can be extended for bimaterial interface
% for center crack use 'BMatrixXIGA_flexCentre'
% can be used for multiple edge crack 

global controlPts element 

sctr = element(e,:);
nn   = length(sctr);

% compute the jacobian of physical and parameter domain mapping
% then the derivative w.r.t spatial physical coordinates

pts        = controlPts(sctr,:);
x_ptss=pts(:,1);
y_ptss=pts(:,2);
jacob      = pts'*[dRdxi' dRdeta'];
J1         = det(jacob);
invJacob   = inv(jacob);
dRdx       = [dRdxi' dRdeta'] * invJacob; % [dRdx dRdy]

Gpt        = N * pts;            % GP in global coord, used   % X and Y as functions of xi and eta
                                   % this our x and y

% natural second derivatives
X_xi=jacob(1,1);
Y_xi=jacob(2,1);
X_eta=jacob(1,2);
Y_eta=jacob(2,2);

bigMatrix =  [X_xi^2 2*X_xi*Y_xi Y_xi^2; X_xi* X_eta (Y_eta*X_xi+Y_xi*X_eta) Y_xi*Y_eta;... % see notes  
                        X_eta^2 2*X_eta*Y_eta Y_eta^2];
INVbigMatrix =inv(bigMatrix);

X_xixi= [dR2dxi dR2dxi;dR2dxideta dR2dxideta;dR2deta dR2deta]*[x_ptss zeros(nn,1);zeros(nn,1) y_ptss];
                                %3by2 matrix in right side
% right side of the equation
dRxx=INVbigMatrix*([dR2dxi;dR2dxideta;dR2deta]-X_xixi*dRdx');  %[Rxx; Rxy; Ryy]

% Bfem is always computed

Bfem = zeros(5,3*nn);
Bfem(1,1:3:3*nn)  = dRdx(:,1)';
Bfem(2,2:3:3*nn)  = dRdx(:,2)';
Bfem(3,1:3:3*nn)  = dRdx(:,2)';
Bfem(3,2:3:3*nn)  = dRdx(:,1)';

Bfem(4,3:3:3*nn)  = dRdx(:,1)';
Bfem(5,3:3:3*nn)  = dRdx(:,2)';



% Switch between non-enriched and enriched elements

if ( any(enrich_node(sctr)) == 0 ) % Non-enriched elements
    B_piezo = Bfem;
else                               % Enriched elements
    Bxfem = [] ;
    % loop on nodes, check node is enriched ...
    for in = 1 : nn
        nodeId = sctr(in);
        enrnoI = enrich_node(nodeId); 
        dRidx  = dRdx(in,1);
        dRidy  = dRdx(in,2);
        Ri     = N(in);
        
        if ( enrnoI == 1)     % H(x) enriched node
            
            
            % Enrichment function, H(x) at global Gauss point
            dist = signed_distance(xCr,Gpt); % signed distance from gauss point (ie is 'x')
            Hgp  = heaviside(dist);
            
            % Enrichment function, H(x) at node "in"
            dist = signed_distance(xCr,controlPts(nodeId,:));
            Hi   = heaviside(dist); % simply signed distance
            %Hi=0;
            % Bxfem at node "in"
            
            aa   = dRidx*(Hgp - Hi);
            bb   = dRidy*(Hgp - Hi);
            
            BI_enr = [aa 0 0;
                      0 bb 0;
                      bb aa 0;
                      0 0 aa;
                      0 0 bb;];
            % Add to the total Bxfem
            Bxfem = [Bxfem BI_enr]; % done for split quadrature ...This will be added for each nodes
            clear BI_enr ;
  %---------------------carck tip %enrichment-----------------------------%------------------------------
  
        elseif ( enrnoI == 2) % B(x) enriched node
            
            seg     = xCr(2,:) - xCr(1,:);   % tip segment
            alpha   = atan2(seg(2),seg(1));  % inclination angle
            QT      = [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)]; % angle transformation
            
            % compute branch functions at Gauss point
            
            xp    = QT*(Gpt-xTip)';           % local coordinates
            r     = sqrt(xp(1)*xp(1)+xp(2)*xp(2));
            theta = atan2(xp(2),xp(1));
            
            if ( theta > pi || theta < -pi)
                disp (['something wrong with angle ',num2str(thet)]);
            end
            
            [Br,dBdx,dBdy] = branch(r,theta,alpha);
            
            % compute branch functions at node "in"
            
            xp    = QT*(controlPts(nodeId,:)-xTip)';
            r     = sqrt(xp(1)*xp(1)+xp(2)*xp(2));
            theta = atan2(xp(2),xp(1));
            
            if ( theta > pi || theta < -pi)
                disp (['something wrong with angle ',num2str(thet)]);
            end
            
            [BrI] = branch_node(r,theta); %BrI = zeros(4);
            
            % components of Benr matrix
            
            aa = dRidx*(Br(1)-BrI(1)) + Ri*dBdx(1) ;
            bb = dRidy*(Br(1)-BrI(1)) + Ri*dBdy(1) ;
            B1_enr = [aa 0 0 ; 0 bb 0 ; bb aa 0;0 0 aa;0 0 bb];
            
            aa = dRidx*(Br(2)-BrI(2)) + Ri*dBdx(2) ;
            bb = dRidy*(Br(2)-BrI(2)) + Ri*dBdy(2) ;
            B2_enr = [aa 0 0 ; 0 bb 0 ; bb aa 0;0 0 aa;0 0 bb];
            
            aa = dRidx*(Br(3)-BrI(3)) + Ri*dBdx(3) ;
            bb = dRidy*(Br(3)-BrI(3)) + Ri*dBdy(3) ;
            B3_enr = [aa 0 0 ; 0 bb 0 ; bb aa 0;0 0 aa;0 0 bb];
            
            aa = dRidx*(Br(4)-BrI(4)) + Ri*dBdx(4) ;
            bb = dRidy*(Br(4)-BrI(4)) + Ri*dBdy(4) ;
            B4_enr = [aa 0 0 ; 0 bb 0 ; bb aa 0;0 0 aa;0 0 bb];
            
            BI_enr = [B1_enr B2_enr B3_enr B4_enr];
            Bxfem = [Bxfem BI_enr];
            clear BI_enr ;
            clear B1_enr; clear B2_enr; clear B3_enr; clear B4_enr;
            
            %--------------%end for edge crack problems-----------
            
            
           
            
            
%-------------------------------------------------------------------------------------------------%------------------
%
        elseif ( enrnoI == 3)     % material interface enriched node
            chi  = CHI(sctr);                        
            Zm   = dot(N,chi);                        
            
            % enrichment functions and derivatives
            
            if     (weakEnrFunc == 1)
                % Moes function
                achi = abs(chi);
                Za   = dot(N,achi);
                Zmn  = Zm/abs(Zm);
                E    = Za-abs(Zm);
                Exi  = dot(dRdxi,achi)  - Zmn*dot(dRdxi,chi);
                Eet  = dot(dRdeta,achi) - Zmn*dot(dRdeta,chi);
            elseif (weakEnrFunc == 2)
                E    = abs(Zm);
                %E    = abs(Zm)-abs(chi(in)); % shift enrichment
                Exi  = sign(Zm)*dot(dRdxi,chi);
                Eet  = sign(Zm)*dot(dRdeta,chi);
            end
            
            aa = dRdxi (in)*E + Ri*Exi;
            bb = dRdeta(in)*E + Ri*Eet;
            
            t  = [aa bb]*invJacob; % derivatives w.r.t global coords
            
            BI_enr = [t(1) 0;0 t(2);t(2) t(1)] ;
            
            % Add to the total Bxfem
            Bxfem = [Bxfem BI_enr];
            clear BI_enr ;
            
        elseif ( enrnoI == 4) % bi-mat crack tip/mat interface enriched node
            crackId = crack_node(nodeId);
            xTip    = xTips(crackId,:);
            xCr     = reshape(xCrack(crackId,:,:),2,2);
            seg     = xCr(2,:) - xCr(1,:);   % tip segment
            alpha   = atan2(seg(2),seg(1));  % inclination angle
            QT      = [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)];
            
            % compute branch functions at Gauss point
            xp    = QT*(Gpt-xTip)';           % local coordinates
            r     = sqrt(xp(1)*xp(1)+xp(2)*xp(2));
            theta = atan2(xp(2),xp(1));
            
            if ( theta > pi || theta < -pi)
                disp (['something wrong with angle ',num2str(thet)]);
            end
            
            [Br,dBdx,dBdy] = biMatCrackBranch(r,theta,ep,alpha);
                      
            % compute branch functions at node "in"
            
            xp    = QT*(controlPts(nodeId,:)-xTip)';
            r     = sqrt(xp(1)*xp(1)+xp(2)*xp(2));
            theta = atan2(xp(2),xp(1));
            
            if ( theta > pi || theta < -pi)
                disp (['something wrong with angle ',num2str(thet)]);
            end
            
            [BrI] = biMatCrackBranchNode(r,theta,ep,alpha);
            
            % components of Benr matrix due to branch functions
            
            for i=1:12
                aa = dRidx*(Br(i)) + Ri*dBdx(i) ;
                bb = dRidy*(Br(i)) + Ri*dBdy(i) ;
                BI_enr = [aa 0 ; 0 bb ; bb aa];
                                
                Bxfem = [Bxfem BI_enr];                                
            end
            
            % enrichment functions and derivatives due to material interface 
            
            chi  = CHI(sctr);                        
            Zm   = dot(N,chi);                        
            
            % enrichment functions and derivatives
            
            if     (weakEnrFunc == 1)
                % Moes function
                achi = abs(chi);
                Za   = dot(N,achi);
                Zmn  = Zm/abs(Zm);
                E    = Za-abs(Zm);
                Exi  = dot(dRdxi,achi)  - Zmn*dot(dRdxi,chi);
                Eet  = dot(dRdeta,achi) - Zmn*dot(dRdeta,chi);
            elseif (weakEnrFunc == 2)
                E    = abs(Zm);
                %E    = abs(Zm)-abs(chi(in)); % shift enrichment
                Exi  = sign(Zm)*dot(dRdxi,chi);
                Eet  = sign(Zm)*dot(dRdeta,chi);
            end
            
            aa = dRdxi (in)*E + Ri*Exi;
            bb = dRdeta(in)*E + Ri*Eet;
            
            t  = [aa bb]*invJacob; 
            
            BI_enr = [t(1) 0;0 t(2);t(2) t(1)] ;
            
            Bxfem = [Bxfem BI_enr]; 
            clear BI_enr ;

         elseif ( enrnoI == 5) % bi-mat crack tip enriched node
            crackId = crack_node(nodeId);
            xTip    = xTips(crackId,:);
            xCr     = reshape(xCrack(crackId,:,:),2,2);
            seg     = xCr(2,:) - xCr(1,:);   % tip segment
            alpha   = atan2(seg(2),seg(1));  % inclination angle
            QT      = [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)];
            
            % compute branch functions at Gauss point
            xp    = QT*(Gpt-xTip)';           % local coordinates
            r     = sqrt(xp(1)*xp(1)+xp(2)*xp(2));
            theta = atan2(xp(2),xp(1));
            
            if ( theta > pi || theta < -pi)
                disp (['something wrong with angle ',num2str(thet)]);
            end
            
            [Br,dBdx,dBdy] = biMatCrackBranch(r,theta,ep,alpha);
                      
            % compute branch functions at node "in"
            
            xp    = QT*(controlPts(nodeId,:)-xTip)';
            r     = sqrt(xp(1)*xp(1)+xp(2)*xp(2));
            theta = atan2(xp(2),xp(1));
            
            if ( theta > pi || theta < -pi)
                disp (['something wrong with angle ',num2str(thet)]);
            end
            
            [BrI] = biMatCrackBranchNode(r,theta,ep,alpha);
            
            % components of Benr matrix due to branch functions
            
            for i=1:12
                aa = dRidx*(Br(i)) + Ri*dBdx(i) ;
                bb = dRidy*(Br(i)) + Ri*dBdy(i) ;
                BI_enr = [aa 0 ; 0 bb ; bb aa];
                                
                Bxfem = [Bxfem BI_enr];                                
            end
        elseif ( enrnoI == 6) % 4 B(x) enriched node and 1 weak enr. func.
            crackId = crack_node(nodeId);
            xTip    = xTips(crackId,:);
            xCr     = reshape(xCrack(crackId,:,:),2,2);
            seg     = xCr(2,:) - xCr(1,:);   % tip segment
            alpha   = atan2(seg(2),seg(1));  % inclination angle
            QT      = [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)];
            
            % compute branch functions at Gauss point
            xp    = QT*(Gpt-xTip)';           % local coordinates
            r     = sqrt(xp(1)*xp(1)+xp(2)*xp(2));
            theta = atan2(xp(2),xp(1));
            
            if ( theta > pi || theta < -pi)
                disp (['something wrong with angle ',num2str(thet)]);
            end
            
            [Br,dBdx,dBdy] = branch(r,theta,alpha);
            
            % compute branch functions at node "in"
            
            xp    = QT*(controlPts(nodeId,:)-xTip)';
            r     = sqrt(xp(1)*xp(1)+xp(2)*xp(2));
            theta = atan2(xp(2),xp(1));
            
            if ( theta > pi || theta < -pi)
                disp (['something wrong with angle ',num2str(thet)]);
            end
            
            [BrI] = branch_node(r,theta); %BrI = zeros(4);
            
            % components of Benr matrix
            
            aa = dRidx*(Br(1)-BrI(1)) + Ri*dBdx(1) ;
            bb = dRidy*(Br(1)-BrI(1)) + Ri*dBdy(1) ;
            B1_enr = [aa 0 ; 0 bb ; bb aa];
            
            aa = dRidx*(Br(2)-BrI(2)) + Ri*dBdx(2) ;
            bb = dRidy*(Br(2)-BrI(2)) + Ri*dBdy(2) ;
            B2_enr = [aa 0 ; 0 bb ; bb aa];
            
            aa = dRidx*(Br(3)-BrI(3)) + Ri*dBdx(3) ;
            bb = dRidy*(Br(3)-BrI(3)) + Ri*dBdy(3) ;
            B3_enr = [aa 0 ; 0 bb ; bb aa];
            
            aa = dRidx*(Br(4)-BrI(4)) + Ri*dBdx(4) ;
            bb = dRidy*(Br(4)-BrI(4)) + Ri*dBdy(4) ;
            B4_enr = [aa 0 ; 0 bb ; bb aa];
            
            BI_enr = [B1_enr B2_enr B3_enr B4_enr];
            Bxfem = [Bxfem BI_enr];
            clear BI_enr ;
            clear B1_enr; clear B2_enr; clear B3_enr; clear B4_enr;
            
            % enrichment functions and derivatives due to material interface 
            
            chi  = CHI(sctr);                        
            Zm   = dot(N,chi);                        
            
            % enrichment functions and derivatives
            
            if     (weakEnrFunc == 1)
                % Moes function
                achi = abs(chi);
                Za   = dot(N,achi);
                Zmn  = Zm/abs(Zm);
                E    = Za-abs(Zm);
                Exi  = dot(dRdxi,achi)  - Zmn*dot(dRdxi,chi);
                Eet  = dot(dRdeta,achi) - Zmn*dot(dRdeta,chi);
            elseif (weakEnrFunc == 2)
                E    = abs(Zm);
                Exi  = sign(Zm)*dot(dRdxi,chi);
                Eet  = sign(Zm)*dot(dRdeta,chi);
            end
            
            aa = dRdxi (in)*E + Ri*Exi;
            bb = dRdeta(in)*E + Ri*Eet;
            
            t  = [aa bb]*invJacob; 
            
            BI_enr = [t(1) 0;0 t(2);t(2) t(1)] ;
            
            Bxfem = [Bxfem BI_enr]; 
            clear BI_enr ;
        end
    end          % end of loop on nodes
    % B matrix
 %%%----------------------------------------&&&&&-------------------------------------------------------&&&&-----------------------------------
    B_piezo = [Bfem Bxfem];
    clear Bfem; clear Bxfem;

    
 
end              % end of switch between enriched and non-enriched elements
%%
%&&&&&&&&&&&&&& B Flexo &&&&&&&&&&&&&&&&&&&&&&&&&
    %B FEM flexo
    
BfemFlex = zeros(8,3*nn); % 6 strain gradients and 2 electric fields

BfemFlexPost=zeros(10,3*nn);  % used only in post processing

% for K matrix

BfemFlex(1,1:3:3*nn)  = dRxx(1,:);
BfemFlex(2,2:3:3*nn)  = dRxx(2,:);
BfemFlex(3,1:3:3*nn)  = dRxx(2,:);
BfemFlex(3,2:3:3*nn)  = dRxx(1,:);

BfemFlex(4,1:3:3*nn)  = dRxx(2,:);
BfemFlex(5,2:3:3*nn)  = dRxx(3,:);
BfemFlex(6,1:3:3*nn)  = dRxx(3,:);
BfemFlex(6,2:3:3*nn)  = dRxx(2,:);

BfemFlex(7,3:3:3*nn)  = dRdx(:,1)';
BfemFlex(8,3:3:3*nn)  = dRdx(:,2)';

% for post processing

BfemFlexPost(1,1:3:3*nn)  = dRxx(1,:);
BfemFlexPost(2,2:3:3*nn)  = dRxx(2,:);
BfemFlexPost(3,1:3:3*nn)  = dRxx(2,:);
BfemFlexPost(3,2:3:3*nn)  = dRxx(1,:);

BfemFlexPost(4,1:3:3*nn)  = dRxx(2,:);
BfemFlexPost(5,2:3:3*nn)  = dRxx(3,:);
BfemFlexPost(6,1:3:3*nn)  = dRxx(3,:);
BfemFlexPost(6,2:3:3*nn)  = dRxx(2,:);

BfemFlexPost(7,3:3:3*nn)  = dRxx(1,:);
BfemFlexPost(8,3:3:3*nn)  = dRxx(2,:);
BfemFlexPost(9,3:3:3*nn)  = dRxx(2,:);
BfemFlexPost(10,3:3:3*nn)  = dRxx(3,:);


% Switch between non-enriched and enriched elements
if ( any(enrich_node(sctr)) == 0 ) % Non-enriched elements
    B_flexo = BfemFlex;
    B_flexoPost =BfemFlexPost;
    
    
    
    % -----------&&&&&& B XFEM flexo &&&&&&-----------------------

else                               % Enriched elements
    BxfemFlexo = [] ;
    BxfemFlexoPost=[];
    
    % loop on nodes, check node is enriched ...
    for in = 1 : nn
        nodeId = sctr(in);
        enrnoI = enrich_node(nodeId);
        dRidx  = dRdx(in,1);
        dRidy  = dRdx(in,2);
        Ri     = N(in);
        dRi2dx=dRxx(1,in);
        dRi2dy=dRxx(3,in);
        dRi2dxdy=dRxx(2,in);
        if ( enrnoI == 1)     % H(x) enriched node
            
           
            
            % Enrichment function, H(x) at global Gauss point
            dist = signed_distance(xCr,Gpt); % signed distance from gauss point (ie is 'x')
            Hgp  = heaviside(dist);
            
            % Enrichment function, H(x) at node "in"
            dist = signed_distance(xCr,controlPts(nodeId,:));
            Hi   = heaviside(dist); % simply signed distance
            %Hi=0;
            % Bxfem at node  ...   "in"
            aa   = dRidx*(Hgp - Hi);
            bb   = dRidy*(Hgp - Hi);
            
            cc=dRi2dx*(Hgp - Hi);
            dd=dRi2dxdy*(Hgp - Hi);
            ee=dRi2dy*(Hgp - Hi);
            
             BI_enrFlexo = [cc 0 0;
                            0 dd 0;
                            dd cc 0;
                            dd 0 0;
                            0 ee 0;
                            ee cc 0;
                            0 0 aa;
                            0 0 bb;];
                        
             %post
              BI_enrFlexopost = [cc 0 0;
                                 0 dd 0;
                                 dd cc 0;
                                 dd 0 0;
                                 0 ee 0;
                                 ee cc 0;
                                 0 0 cc;
                                 0 0 dd;
                                 0 0 dd;
                                 0 0 ee];
                        
            % Add to the total BxfemFlexo

           BxfemFlexo = [BxfemFlexo BI_enrFlexo]; % done for split nodes ...This will be added for each nodes
           BxfemFlexoPost=[BxfemFlexoPost BI_enrFlexopost];
            clear BI_enrFlexo ;
           clear BI_enrFlexopost;
            
            %----------Crack Tip enrichment for flexoelectric materials-------------%
            
            elseif ( enrnoI == 2) % B(x) enriched node
            
            seg     = xCr(2,:) - xCr(1,:);   % tip segment
            alpha   = atan2(seg(2),seg(1));  % inclination angle
            QT      = [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)]; % angle transformation
            
            % compute branch functions at Gauss point
            
            xp    = QT*(Gpt-xTip)';           % local coordinates
            r     = sqrt(xp(1)*xp(1)+xp(2)*xp(2));
            theta = atan2(xp(2),xp(1));
            
            if ( theta > pi || theta < -pi)
                disp (['something wrong with angle ',num2str(thet)]);
            end
            
            [Br,dBdx,dBdy,dBdxx,dBdyy,dBdxy] = branch(r,theta,alpha);
            
             % compute branch functions at node "in"
            
            xp    = QT*(controlPts(nodeId,:)-xTip)';
            r     = sqrt(xp(1)*xp(1)+xp(2)*xp(2));
            theta = atan2(xp(2),xp(1));
            
            if ( theta > pi || theta < -pi)
                disp (['something wrong with angle ',num2str(thet)]);
            end
            
            [BrI] = branch_node(r,theta); %BrI = zeros(4);
            
                      % components of Benr matrix
            
            aa = dRidx*(Br(1)-BrI(1)) + Ri*dBdx(1) ;
            bb = dRidy*(Br(1)-BrI(1)) + Ri*dBdy(1) ;
            
            cc = dRi2dx*(Br(1)-BrI(1)) + 2*dBdx(1)*dRidx + Ri*dBdxx(1);  % xx
            dd = dRi2dy*(Br(1)-BrI(1)) + 2*dBdy(1)*dRidy + Ri*dBdyy(1);  %yy
            ee = dRi2dxdy*(Br(1)-BrI(1)) + dBdy(1)*dRidx +dBdx(1)*dRidy + Ri*dBdxy(1);  %xy
            
            B1_enr = [cc 0 0 ; 0 ee 0 ; ee cc 0;ee 0 0;0 dd 0;dd ee 0;0 0 aa;0 0 bb];
            B1_enrPost = [cc 0 0 ; 0 ee 0 ; ee cc 0;ee 0 0;0 dd 0;dd ee 0;0 0 cc;0 0 ee;0 0 ee;0 0 dd];
            %-----
            aa = dRidx*(Br(2)-BrI(2)) + Ri*dBdx(2) ;
            bb = dRidy*(Br(2)-BrI(2)) + Ri*dBdy(2) ;
            
            cc = dRi2dx*(Br(2)-BrI(2)) + 2*dBdx(2)*dRidx + Ri*dBdxx(2);  % xx
            dd = dRi2dy*(Br(2)-BrI(2)) + 2*dBdy(2)*dRidy + Ri*dBdyy(2);  %yy
            ee = dRi2dxdy*(Br(2)-BrI(2)) + dBdy(2)*dRidx +dBdx(2)*dRidy + Ri*dBdxy(2);  %xy
            
            B2_enr = [cc 0 0 ; 0 ee 0 ; ee cc 0;ee 0 0;0 dd 0;dd ee 0;0 0 aa;0 0 bb];
            B2_enrPost = [cc 0 0 ; 0 ee 0 ; ee cc 0;ee 0 0;0 dd 0;dd ee 0;0 0 cc;0 0 ee;0 0 ee;0 0 dd];
            %-----
             aa = dRidx*(Br(3)-BrI(3)) + Ri*dBdx(3) ;
            bb = dRidy*(Br(3)-BrI(3)) + Ri*dBdy(3) ;
            
            cc = dRi2dx*(Br(3)-BrI(3)) + 2*dBdx(3)*dRidx + Ri*dBdxx(3);  % xx
            dd = dRi2dy*(Br(3)-BrI(3)) + 2*dBdy(3)*dRidy + Ri*dBdyy(3);  %yy
            ee = dRi2dxdy*(Br(3)-BrI(3)) + dBdy(3)*dRidx +dBdx(3)*dRidy + Ri*dBdxy(3);  %xy
            
            B3_enr = [cc 0 0 ; 0 ee 0 ; ee cc 0;ee 0 0;0 dd 0;dd ee 0;0 0 aa;0 0 bb];
            B3_enrPost = [cc 0 0 ; 0 ee 0 ; ee cc 0;ee 0 0;0 dd 0;dd ee 0;0 0 cc;0 0 ee;0 0 ee;0 0 dd];
            %---------
             aa = dRidx*(Br(4)-BrI(4)) + Ri*dBdx(4) ;
            bb = dRidy*(Br(4)-BrI(4)) + Ri*dBdy(4) ;
            
            cc = dRi2dx*(Br(4)-BrI(4)) + 2*dBdx(4)*dRidx + Ri*dBdxx(4);  % xx
            dd = dRi2dy*(Br(4)-BrI(4)) + 2*dBdy(4)*dRidy + Ri*dBdyy(4);  %yy
            ee = dRi2dxdy*(Br(4)-BrI(4)) + dBdy(4)*dRidx +dBdx(4)*dRidy + Ri*dBdxy(4);  %xy
            
            B4_enr = [cc 0 0 ; 0 ee 0 ; ee cc 0;ee 0 0;0 dd 0;dd ee 0;0 0 aa;0 0 bb];
            B4_enrPost = [cc 0 0 ; 0 ee 0 ; ee cc 0;ee 0 0;0 dd 0;dd ee 0;0 0 cc;0 0 ee;0 0 ee;0 0 dd];
            
            BI_enrFlexo = [B1_enr B2_enr B3_enr B4_enr];
            BI_enrFlexopost = [B1_enrPost B2_enrPost B3_enrPost B4_enrPost];
            
            BxfemFlexo = [BxfemFlexo BI_enrFlexo];
            BxfemFlexoPost=[BxfemFlexoPost BI_enrFlexopost];
            clear BI_enrFlexo ;
            clear B1_enr; clear B2_enr; clear B3_enr; clear B4_enr;
             clear BI_enrFlexopost ;
            clear B1_enrPost; clear B2_enrPost; clear B3_enrPost; clear B4_enrPost;
        end
        end
    B_flexo = [BfemFlex BxfemFlexo];
    B_flexoPost = [BfemFlexPost BxfemFlexoPost];
    clear BfemFlex; clear BxfemFlexo;
    clear BfemFlexPost; clear BxfemFlexoPost;
    
end
    

