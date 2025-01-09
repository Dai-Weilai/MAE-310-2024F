clear all; clc;

kappa = 1.0; % conductivity
P_ratio=0.3;
E_modulus=1E9;
Problem_type=input('please enter the problem type! 1 for plane stress and 2 for plane strain!');

% exact solution
exact = @(x,y) x*(1-x)*y*(1-y);
exact_x = @(x,y) (1-2*x)*y*(1-y);
exact_y = @(x,y) x*(1-x)*(1-2*y);

f = @(x,y) 2.0*kappa*x*(1-x) + 2.0*kappa*y*(1-y); % source term

% quadrature rule
n_int_xi  = 3;
n_int_eta = 3;
n_int     = n_int_xi * n_int_eta;
[xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);

% read the mesh 
n_en   = 4;               % number of nodes in an element
coor=zeros(2,n_en);
[n_el,n_np,coor,IEN]=readMshFile('quarter-plate-with-hole-quad.msh');
% generate the nodal coordinates
% for ny = 1 : n_np_y
%   for nx = 1 : n_np_x
%     index = (ny-1)*n_np_x + nx; % nodal index
%     coor(1,index) = (nx-1) * hx;
%     coor(2,index) = (ny-1) * hy;
%   end
% end

% IEN array
% IEN = zeros(n_el, n_en);
% for ex = 1 : n_el_x
%   for ey = 1 : n_el_y
%     ee = (ey-1) * n_el_x + ex; % element index
%     IEN(ee, 1) = (ey-1) * n_np_x + ex;
%     IEN(ee, 2) = (ey-1) * n_np_x + ex + 1;
%     IEN(ee, 3) =  ey    * n_np_x + ex + 1;
%     IEN(ee, 4) =  ey    * n_np_x + ex;
%   end
% end

% ID array
ID = zeros(n_np,2);
counter = 0;
for i=1:n_np
    if coor(i,1) == -1 
        ID(1,i)=65535;
        elseif coor(i,2) == -1
        ID(2,i)=65535;    
    end
end
count=0;
for i=1:n_np
    for j=1:2
        if ID(j,i) == 0
            count =count+1;
            ID(j,i)=count;
        else 
            ID(j,i)=0;
        end
    end
end
n_eq=count;
% for ny = 2 : n_np_y - 1
%   for nx = 2 : n_np_x - 1
%     index = (ny-1)*n_np_x + nx;
%     counter = counter + 1;
%     ID(index) = counter;  
%   end
% end
% 
% n_eq = counter;

LM = ID(IEN);

% allocate the stiffness matrix and load vector
K = spalloc(n_eq, n_eq, 9 * n_eq);
F = zeros(n_eq, 1);


% loop over element to assembly the matrix and vector
for ee = 1 : n_el
  x_ele = x_coor( IEN(ee, 1:n_en) );
  y_ele = y_coor( IEN(ee, 1:n_en) );
  
  k_ele = zeros(n_en, n_en); % element stiffness matrix
  f_ele = zeros(n_en, 1);    % element load vector
  
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
      
      f_ele(aa) = f_ele(aa) + weight(ll) * detJ * f(x_l, y_l) * Na;
      
      for bb = 1 : n_en
        Nb = Quad(bb, xi(ll), eta(ll));
        [Nb_xi, Nb_eta] = Quad_grad(bb, xi(ll), eta(ll));
        Nb_x = (Nb_xi * dy_deta - Nb_eta * dy_dxi) / detJ;
        Nb_y = (-Nb_xi * dx_deta + Nb_eta * dx_dxi) / detJ;
        
        k_ele(aa, bb) = k_ele(aa,bb) + weight(ll) * detJ * kappa * (Na_x * Nb_x + Na_y * Nb_y);
      end % end of bb loop
    end % end of aa loop
  end % end of quadrature loop
 
  for aa = 1 : n_en
    PP = LM(ee, aa);
    if PP > 0
      F(PP) = F(PP) + f_ele(aa);
      
      for bb = 1 : n_en
        QQ = LM(ee, bb);
        if QQ > 0
          K(PP, QQ) = K(PP, QQ) + k_ele(aa, bb);
        else
          % modify F with the boundary data
          % here we do nothing because the boundary data g is zero or
          % homogeneous
        end
      end  
    end
  end
end

% solve the stiffness matrix
dn = K \ F;

% insert dn back into the vector for all nodes
disp = zeros(n_np, 1);

for ii = 1 : n_np
  index = ID(ii);
  if index > 0
    disp(ii) = dn(index);
  else
    % modify disp with the g data. Here it does nothing because g is zero
  end
end

% save the solution vector and number of elements to disp with name
% HEAT.mat
save("HEAT", "disp", "n_el_x", "n_el_y");

% EOF