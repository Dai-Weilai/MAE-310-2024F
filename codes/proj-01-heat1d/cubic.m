% clear all; clc; clf;
u_exact=@(x) x.^5;
ux_exact=@(x) 5*x.^4;
f = @(x) -20*x.^3; % f(x) is the source
g=1;h=0;% u=g at x=1;-u; x=h at x=0
pp=3; % polynomial degree
n_en = pp + 1;         % number of element or local nodes
n_int=10;
% Setup the quadrature rule
[xi, weight] = Gauss(n_int, -1, 1);
error_resultL=zeros(1,8);
error_resultH=zeros(1,8);
ref=zeros(1,8);
Lfenzi=zeros(8);
Lfenmu=zeros(8);
Hfenzi=zeros(8);
Hfenmu=zeros(8);
for n_el=2:2:16 % number of elements
  % Setup the mesh
  n_np = n_el * pp + 1;  % number of nodal points
  n_eq = n_np - 1;       % number of equations  
  hh=1/n_eq;
  x_coor = 0 : hh : 1;   % nodal coordinates for equally spaced nodes

  % Setup the IEN
    IEN = zeros(n_el, n_en);
    for ee = 1 : n_el
        for aa = 1 : n_en
            IEN(ee, aa) = (ee - 1) * pp + aa;
        end
    end
  % Setup the ID array
    ID = 1 : n_np;
    ID(end) = 0;
  % allocate the stiffness matrix
  K = spalloc(n_eq, n_eq, (2*pp+1)*n_eq);
  F = zeros(n_eq, 1);
  % Assembly of the stiffness matrix and load vector
  for ee = 1 : n_el
    k_ele = zeros(n_en, n_en); % allocate a zero element stiffness matrix
    f_ele = zeros(n_en, 1);    % allocate a zero element load vector
    x_ele = x_coor(IEN(ee,:)); % x_ele(aa) = x_coor(A) with A = IEN(aa, ee)
    % quadrature loop
    for qua = 1 : n_int    
      dx_dxi = 0.0;
      x_l = 0.0;
      for aa = 1 : n_en
        x_l    = x_l    + x_ele(aa) * PolyShape(pp, aa, xi(qua), 0);
        dx_dxi = dx_dxi + x_ele(aa) * PolyShape(pp, aa, xi(qua), 1);
      end
      dxi_dx = 1.0 / dx_dxi;
      for aa = 1 : n_en
        f_ele(aa) = f_ele(aa) + weight(qua) * PolyShape(pp, aa, xi(qua), 0) * f(x_l) * dx_dxi;
        for bb = 1 : n_en
          k_ele(aa, bb) = k_ele(aa, bb) + weight(qua) * PolyShape(pp, aa, xi(qua), 1) * PolyShape(pp, bb, xi(qua), 1) * dxi_dx;
        end
      end
    end
    % Assembly of the matrix and vector based on the ID or LM data
    for aa = 1 : n_en
      P = ID(IEN(ee,aa));
      if(P > 0)
        F(P) = F(P) + f_ele(aa);
        for bb = 1 : n_en
          Q = ID(IEN(ee,bb));
          if(Q > 0)
            K(P, Q) = K(P, Q) + k_ele(aa, bb);
          else
            F(P) = F(P) - k_ele(aa, bb) * g; % handles the Dirichlet boundary data
          end
        end
      end
    end
  end
    % ee = 1 F = NA(0)xh
    F(ID(IEN(1,1))) = F(ID(IEN(1,1))) + h;
    % Solve Kd = F equation
    d_temp = K \ F;
    disp = [d_temp; g];
    n_sam = 20;
    xi_sam = -1 : (2/n_sam) : 1;
    x_sam = zeros(n_el * n_sam + 1, 1);
    y_sam = x_sam; % store the exact solution value at sampling points
    u_sam = x_sam; % store the numerical solution value at sampling pts

    for ee = 1 : n_el
      x_ele = x_coor( IEN(ee, :) );
      u_ele = disp( IEN(ee, :) );
      if ee == n_el
        n_sam_end = n_sam+1;
      else
        n_sam_end = n_sam;
      end

      for ll = 1 : n_sam_end
        x_l = 0.0;
        u_l = 0.0;
        for aa = 1 : n_en
          x_l = x_l + x_ele(aa) * PolyShape(pp, aa, xi_sam(ll), 0);
          u_l = u_l + u_ele(aa) * PolyShape(pp, aa, xi_sam(ll), 0);
        end

         x_sam( (ee-1)*n_sam + ll ) = x_l;
         u_sam( (ee-1)*n_sam + ll ) = u_l;
         y_sam( (ee-1)*n_sam + ll ) = x_l^5;
      end
      for i=1:n_int
          x_l=0;
          uh=0;
          for j=1:n_en
              x_l=x_l+x_ele(j)*PolyShape(pp,j,xi(i),0);
              uh=uh+u_ele(j)*PolyShape(pp,j,xi(i),0);
          end 
          uhx=0;
          dx_dxi=0;
          for j=1:n_en
            dx_dxi = dx_dxi + x_ele(j) * PolyShape(pp, j, xi(i), 1);
            uhx=uhx+u_ele(j)* PolyShape(pp, j, xi(i), 1);
          end
        Lfenmu(ee)=Lfenmu(ee)+(u_exact(x_l)^2*dx_dxi)*weight(i);
        Lfenzi(ee)=Lfenzi(ee)+((uh-u_exact(x_l))^2*dx_dxi)*weight(i);
        Hfenmu(ee)=Hfenmu(ee)+(ux_exact(x_l)^2*dx_dxi)*weight(i);
        Hfenzi(ee)=Hfenzi(ee)+((uhx-ux_exact(x_l))^2*dx_dxi)*weight(i);
      end
    end
sumL=sum(Lfenzi)/sum(Lfenmu);
sumH=sum(Hfenzi)/sum(Hfenmu);
error_resultL(n_el/2)=sqrt(sumL);
error_resultH(n_el/2)=sqrt(sumH);
ref(n_el/2)=hh;
end
x=[2,4,6,8,10,12,14,16];
plot(log10(ref),log10(error_resultH));
hold on;
plot(log10(ref),log10(error_resultL));