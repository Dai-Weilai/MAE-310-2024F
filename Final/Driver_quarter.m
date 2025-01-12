clear all; 
clc;
P_ratio=0.3;
E_modulus=10^9;
Shear_m=E_modulus/(2+2*P_ratio);
Lam = (E_modulus*P_ratio)/((1+P_ratio)*(1-2*P_ratio));
% Problem_type=input('please enter the problem type! 1 for plane stress and 2 for plane strain!');
% if Problem_type ==1
%     D=(E_modulus/(1-P_ratio^2))*[1,P_ratio,0;P_ratio,1,0;0,0,(1 - P_ratio)/2];%stress
% elseif Problem_type ==2
    D=[2*Shear_m+Lam,Lam,0;Lam,2*Shear_m+Lam,0;0,0,Shear_m];%strain
% end

addpath("mesh\");
% read the mesh 
n_en   = 4;               % number of nodes in an element
coor=zeros(500,2);
[n_np,n_el,coor,IEN]=readMshFile('quarter-plate-with-hole-quad.msh');%其他refine过1至3遍的msh文件也在目录中,可以自行取用！

%exact solution
sigmarr=@(Tx,thet,R,r) Tx./2.*(1 - R.^2./r.^2)+Tx./2.*(1 - 4.*R.^2./r.^2 + 3.*R.^4./r.^4).*cos(2.*thet);
sigmart=@(Tx,thet,R,r) -Tx./2.*(1 + 2.*R.^2./r.^2 - 3.*R.^4./r.^4).*sin(2.*thet);
sigmatt=@(Tx,thet,R,r) Tx./2.*(1 + R.^2./r.^2) - Tx./2.*(1 + 3.*R.^4./r.^4).*cos(2.*thet);

Tx=10000;
R=0.3;

% Boundry condition
% Dirichle=input('Where to apply the BC? Please input in the order of: left,right,Above,Below.0 for not apply,1 for apply.');%你可以自行选择在哪几条边上施加边界。
Dirichle=[1,0,0,1];
BC=zeros(n_np,2);
if Dirichle(1)==1
    for i=1:n_np
        if coor(i,1)==-1
            BC(i)=coor(i);
        end
    end
end
if Dirichle(2)==1
    for i=1:n_np
        if coor(i,1)==1
            BC(i)=coor(i);
        end
    end
end
if Dirichle(3)==1
    for i=1:n_np
        if coor(i,2)==1
            BC(i)=coor(i);
        end
    end
end
if Dirichle(4)==1
    for i=1:n_np
        if coor(i,2)==-1
            BC(i)=coor(i);
        end
    end
end
r=sqrt((BC(:,1)-(-1)).^2+(BC(:,2)-(-1)).^2);
thet=atan((BC(:,2)-(-1))./((BC(:,1)-(-1))));
sigma_rr=sigmarr(Tx,thet,R,r);
sigma_tt=sigmatt(Tx,thet,R,r);
sigma_rt=sigmart(Tx,thet,R,r);

g=0;
% g=input('please enter the Neumann BC g!');

sigma_m=(sigma_rr+sigma_tt).*0.5;
sigma_Rad=sqrt(((sigma_rr+sigma_tt).*0.5).^2+sigma_rt.^2);

fx=sigma_m+sigma_Rad;
fy=sigma_m-sigma_Rad;
fxy=sigma_Rad;


% ID array
ID = zeros(n_np,2);
counter = 0;
for i=1:n_np
    if ((coor(i,1)==-1)&&(Dirichle(1)==1))||((coor(i,1)==1)&&(Dirichle(2)==1))
        ID(1,i)=65535;
        elseif ((coor(i,2)==-1)&&(Dirichle(3)==1))||((coor(i,2)==1)&&(Dirichle(4)==1))
        ID(2,i)=65535;    
    end %标记ID array中为0的点
end
count=0;
for i=1:n_np
    for j=1:2
        if ID(i,j) == 0
            count =count+1;
            ID(i,j)=count;
        else 
            ID(i,j)=0;
        end
    end
end
n_eq=count;

%LM array
LM =zeros(n_el,n_en,2);
ID_row_1=ID(:,1); ID_row_2=ID(:,2); 
LM(:,:,1)=ID_row_1(IEN);
LM(:,:,2)=ID_row_2(IEN);%两个自由度——LM array需要分成两份。


% quadrature rule
n_int_xi  = 3;
n_int_eta = 3;
n_int     = n_int_xi * n_int_eta;
[xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);


% allocate the stiffness matrix and load vector
K = spalloc(n_eq, n_eq, 9 * n_eq);
F = zeros(n_eq, 1);


% loop over element to assembly the matrix and vector
for ee = 1 : n_el
  x_ele = coor( IEN(ee, 1:n_en),1 );
  y_ele = coor( IEN(ee, 1:n_en),2 );
  
  k_ele = zeros(n_en, n_en,2,2); % element stiffness matrix,每个元素因为有2个自由度的缘故被分割为2x2
  f_ele = zeros(n_en, 2);    % element load vector
  
  for ll = 1 : n_int
    x_l = 0.0; y_l = 0.0;
    dx_dxi = 0.0; dx_deta = 0.0;
    dy_dxi = 0.0; dy_deta = 0.0;
    for aa = 1 : n_en
      x_l = x_l + x_ele(aa) * Quad(aa, xi(ll), eta(ll));
      y_l = y_l + y_ele(aa) * Quad(aa, xi(ll), eta(ll));    
      [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
      dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
      dx_deta = dx_deta + x_ele(aa) * Na_eta;
      dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
      dy_deta = dy_deta + y_ele(aa) * Na_eta;
    end
    
    detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;
    
    for aa = 1 : n_en
      Na = Quad(aa, xi(ll), eta(ll));
      [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
      Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
      Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;
      
      f_ele(aa,1)=f_ele(aa,1)+Na*weight(ll)*detJ*fx(IEN(ee,aa));
      f_ele(aa,2)=f_ele(aa,2)+Na*weight(ll)*detJ*fy(IEN(ee,aa));

      B=zeros(3,2);%获取B矩阵
      for bb = 1 : n_en
        Nb = Quad(bb, xi(ll), eta(ll));
        [Nb_xi, Nb_eta] = Quad_grad(bb, xi(ll), eta(ll));
        Nb_x = (Nb_xi * dy_deta - Nb_eta * dy_dxi) / detJ;
        Nb_y = (-Nb_xi * dx_deta + Nb_eta * dx_dxi) / detJ;
        B(1,2)=Nb_x;B(2,2)=Nb_y;B(3,2)=Nb_x;B(3,1)=Nb_y;
        k_ele(aa,bb,:,:)=B'*D*B;
        % k_ele(aa, bb) = k_ele(aa,bb) + weight(ll) * detJ * kappa * (Na_x * Nb_x + Na_y * Nb_y);
      end % end of bb loop
    end % end of aa loop
  end % end of quadrature loop
 
  for aa = 1 : n_en
    PP_1=LM(ee,aa,1);
    PP_2=LM(ee,aa,2);
    if PP_1> 0
      F(PP_1)=F(PP_1)+f_ele(aa,1);
      for bb = 1 : n_en
        QQ_1= LM(ee, bb,1);
         QQ_2= LM(ee, bb,2);
        if QQ_1 > 0
          K(PP_1, QQ_1) = K(PP_1, QQ_1) + k_ele(aa,bb,1,1);
        elseif QQ_2>0
              K(PP_1, QQ_2) = K(PP_1, QQ_2) + k_ele(aa,bb,1,2);  
          % modify F with the boundary data
          % disp = [d_temp; g];
        end
      end  
    end
    if PP_2> 0
      F(PP_2)=F(PP_2)+f_ele(aa,2);
      for bb = 1 : n_en
        QQ_1= LM(ee, bb,1);
         QQ_2= LM(ee, bb,2);
        if QQ_1 > 0
          K(PP_2, QQ_1) = K(PP_2, QQ_1) + k_ele(aa,bb,2,1);
        elseif QQ_2>0
              K(PP_2, QQ_2) = K(PP_2, QQ_2) + k_ele(aa,bb,2,2);  
          % modify F with the boundary data
          % disp = [d_temp; g];
        end
      end  
    end
  end
end

% solve the stiffness matrix
dn = K \ F;

% insert dn back into the vector for all nodes
disp = zeros(n_np, 2);

for ii = 1 : n_np
    for jj=1:2
     index = ID(ii,jj);
     if index > 0
       disp(ii,jj) = dn(index);
     else
       % disp = [d_temp; g]; nothing because g=0!
     end
    end
end



% EOF