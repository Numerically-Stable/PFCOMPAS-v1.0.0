function [N,dNdxi]=lagrange_basis(type,coord,dim)
% ============================================================
% lagrange_basis - Evaluates shape functions and local derivatives
%
% PHASE-FIELD PHYSICS RELEVANCE:
% The exact analytical solution for a 1D phase-field crack profile 
% follows an exponential decay: phi(x) = exp(-|x|/l0). 
% Standard linear elements (like L2, T3, Q4) only provide $C^0$ continuous, 
% linear interpolation. To capture the exponential curve of the damage 
% variable without severe discretization error, linear elements require 
% an extremely dense, heavily refined mesh inside the crack band.
%
% By providing a robust library of higher-order elements (e.g., quadratic 
% Q8/Q9, cubic Q16, or higher-order triangles T6/T10), researchers can 
% capture the steep phase-field gradients much more accurately with fewer 
% degrees of freedom. This also dramatically improves the accuracy of the 
% Verfürth a posteriori error estimators, which rely on resolving 
% higher-order flux jumps.
%
% INPUT:
%   type  : Topological class (e.g., 'Q4', 'T6', 'H20')
%   coord : Parent/Reference coordinates (xi, eta, zeta)
%   dim   : Spatial dimensions
%
% OUTPUT:
%   N     : Shape function vector evaluated at (xi, eta, zeta)
%   dNdxi : Matrix of shape function derivatives w.r.t parent coordinates
%
% Originally written by Jack Chessa (Northwestern University)
% Augmented with high-order elements by Pranjal Saxena (IIT Kanpur)
% ============================================================
    
  if ( nargin == 2 )
    dim=1;
  end
  
  switch type
   
   %% ==========================================================
   % 1D LINE ELEMENTS (Reference domain: [-1, 1])
   % ==========================================================
   case 'L2'  
    % 1---------2
     xi=coord(1);
     N=[(1-xi)/2;(1+xi)/2];
     dNdxi=[-1/2;1/2];
     
   case 'L3' 
    % 1---------2----------3
     xi=coord(1);
     N=[(1-xi)*xi/(-2);1-xi^2;(1+xi)*xi/2];
     dNdxi=[xi-.5;-2*xi;xi+.5];
     
   case 'L4' 
    % 1-----2-----3-----4
     xi=coord(1);
     N=[-9*(xi+1/3)*(xi-1/3)*(xi-1)/16;
         27*(xi+1)*(xi-1/3)*(xi-1)/16;
         -27*(xi+1)*(xi+1/3)*(xi-1)/16;
        9*(xi+1/3)*(xi-1/3)*(xi+1)/16];
     dNdxi=[- (9*(xi - 1)*(xi - 1/3))/16 - ((9*xi + 3)*(xi - 1))/16 - ((9*xi + 3)*(xi - 1/3))/16;
            (27*(xi - 1)*(xi - 1/3))/16 + ((27*xi + 27)*(xi - 1))/16 + ((27*xi + 27)*(xi - 1/3))/16;
            - (27*(xi - 1)*(xi + 1/3))/16 - ((27*xi + 27)*(xi - 1))/16 - ((27*xi + 27)*(xi + 1/3))/16;
            (9*(xi + 1)*(xi - 1/3))/16 + ((9*xi + 3)*(xi + 1))/16 + ((9*xi + 3)*(xi - 1/3))/16];
            
   %% ==========================================================
   % 2D TRIANGULAR ELEMENTS (Reference domain: Area coordinates)
   % ==========================================================
   case 'T3'
    %    1--------------------2
    if size(coord,2) < 2
      disp('Error two coordinates needed for the T3 element')
    else
      xi=coord(1); eta=coord(2);
      N=[1-xi-eta;xi;eta];
      dNdxi=[-1,-1;1,0;0,1];
    end
    
   case 'T3fs'
    if size(coord,2) < 2
      disp('Error two coordinates needed for the T3fs element')
    else
      xi=coord(1); eta=coord(2);
      N=[1-xi-eta;xi;eta];
      dNdxi=[-1,-1;1,0;0,1];
    end
        
   case 'T4'
    if size(coord,2) < 2
      disp('Error two coordinates needed for the T4 element')
    else
      xi=coord(1); eta=coord(2);
      N=[1-xi-eta-3*xi*eta;xi*(1-3*eta);eta*(1-3*xi);9*xi*eta];
      dNdxi=[-1-3*eta,-1-3*xi;
	     1-3*eta, -3*xi;
	     -3*eta,   1-3*xi;
	     9*eta,   9*xi ];
    end
    
   case 'T6'
    % 6-Node Quadratic Triangle
    if size(coord,2) < 2
      disp('Error two coordinates needed for the T6 element')
    else
      xi=coord(1); eta=coord(2);
      N=[1-3*(xi+eta)+4*xi*eta+2*(xi^2+eta^2);
                                  xi*(2*xi-1);
                                eta*(2*eta-1);
                              4*xi*(1-xi-eta);
                                     4*xi*eta;
                              4*eta*(1-xi-eta)];
        
      dNdxi=[4*(xi+eta)-3   4*(xi+eta)-3;
                   4*xi-1              0; 
                        0        4*eta-1;
           4*(1-eta-2*xi)          -4*xi;
                    4*eta           4*xi;
                   -4*eta  4*(1-xi-2*eta)];
    end
    
    case 'T10'
    % 10-Node Cubic Triangle (Added by Pranjal Saxena)
    if size(coord,2) < 2
      disp('Error two coordinates needed for the T10 element')
    else
      xi=coord(1); eta=coord(2);
      N=[(-9/2)*(xi+eta-1/3)*(xi+eta-2/3)*(xi+eta-1);
            (9/2)*xi*(xi-1/3)*(xi-2/3);
            (9/2)*eta*(eta-1/3)*(eta-2/3);
            (27/2)*xi*(xi+eta-2/3)*(xi+eta-1);
            (-27/2)*xi*(xi-1/3)*(xi+eta-1);
            (27/2)*xi*eta*(xi-1/3);
            (27/2)*xi*eta*(eta-1/3);
            (-27/2)*eta*(eta-1/3)*(xi+eta-1);
            (27/2)*eta*(xi+eta-2/3)*(xi+eta-1);
            (-27)*xi*eta*(xi+eta-1)];
        
      dNdxi=[-(27*eta^2)/2-27*eta*xi+18*eta-(27*xi^2)/2+18*xi-(11/2), -(27*eta^2)/2-27*eta*xi+18*eta-(27*xi^2)/2+18*xi-(11/2);
          (27*xi^2)/2 - 9*xi + 1                                    , 0;
          0                                                         , (27*eta^2)/2-9*eta+1;
          (27*eta^2)/2+54*eta*xi-(45*eta)/2+(81*xi^2)/2-45*xi+9     , (9*xi*(6*eta + 6*xi - 5))/2;
          (9*eta)/2+36*xi-27*eta*xi-(81*xi^2)/2-(9/2)               , -(27*xi*(xi-1/3))/2;
          (9*eta*(6*xi-1))/2                                        ,  (27*xi*(xi-1/3))/2;
           (27*eta*(eta-1/3))/2                                     ,  (9*xi*(6*eta-1))/2;
          -(27*eta*(eta-1/3))/2                                     , 36*eta+(9*xi)/2-27*eta*xi-(81*eta^2)/2-(9/2);
          (9*eta*(6*eta+6*xi-5))/2                                  , (81*eta^2)/2+54*eta*xi-45*eta+(27*xi^2)/2-(45*xi)/2+9;
          -27*eta*(eta+2*xi-1)                                      , -27*xi*(2*eta+xi-1)];
    end   
    
   %% ==========================================================
   % 2D QUADRILATERAL ELEMENTS (Reference domain: [-1, 1] x [-1, 1])
   % ==========================================================
   case 'Q4'
    % 4-Node Bilinear Quadrilateral
    if size(coord,2) < 2
      disp('Error two coordinates needed for the Q4 element')
    else
      xi=coord(1); eta=coord(2);
      N=1/4*[ (1-xi)*(1-eta);
              (1+xi)*(1-eta);
              (1+xi)*(1+eta);
              (1-xi)*(1+eta)];
      dNdxi=1/4*[-(1-eta), -(1-xi);
		         1-eta,    -(1+xi);
		         1+eta,      1+xi;
                -(1+eta),   1-xi];
    end
    
    case 'Q8'
    % 8-Node Serendipity Quadrilateral
    if size(coord,2) < 2
      disp('Error two coordinates needed for the Q8 element')
    else
        xi=coord(1); eta=coord(2);
        N=  [-0.25*(1-xi)*(1-eta)*(1+xi+eta);
             -0.25*(1+xi)*(1-eta)*(1-xi+eta);
             -0.25*(1+xi)*(1+eta)*(1-xi-eta);
             -0.25*(1-xi)*(1+eta)*(1+xi-eta);
              0.5*(1-xi^2)*(1-eta);
              0.5*(1+xi)*(1-eta^2);
              0.5*(1-xi^2)*(1+eta);
              0.5*(1-xi)*(1-eta^2)];
        
 dNdxi=[0.25*(1-eta)*(2*xi+eta), 0.25*(1-xi)*(xi+2*eta);
        0.25*(1-eta)*(2*xi-eta), -0.25*(1+xi)*(xi-2*eta);
        0.25*(1+eta)*(2*xi+eta), 0.25*(1+xi)*(xi+2*eta);
        0.25*(1+eta)*(2*xi-eta), -0.25*(1-xi)*(xi-2*eta);
        -(1-eta)*xi, -0.5*(1-xi^2);
        0.5*(1-eta^2), -(1+xi)*eta;
        -(1+eta)*xi, 0.5*(1-xi^2);
        -0.5*(1-eta^2), -(1-xi)*eta]; 
    end
    
   case 'Q9'
    % 9-Node Biquadratic Quadrilateral
    if size(coord,2) < 2
      disp('Error two coordinates needed for the Q9 element')
    else
      xi=coord(1); eta=coord(2);
      N=1/4*[xi*eta*(xi-1)*(eta-1);
             xi*eta*(xi+1)*(eta-1);
             xi*eta*(xi+1)*(eta+1);
             xi*eta*(xi-1)*(eta+1);
            -2*eta*(xi+1)*(xi-1)*(eta-1);
            -2*xi*(xi+1)*(eta+1)*(eta-1);
            -2*eta*(xi+1)*(xi-1)*(eta+1);
            -2*xi*(xi-1)*(eta+1)*(eta-1);
             4*(xi+1)*(xi-1)*(eta+1)*(eta-1)];
      dNdxi=1/4*[eta*(2*xi-1)*(eta-1),xi*(xi-1)*(2*eta-1);
                 eta*(2*xi+1)*(eta-1),xi*(xi+1)*(2*eta-1);
                 eta*(2*xi+1)*(eta+1),xi*(xi+1)*(2*eta+1);
                 eta*(2*xi-1)*(eta+1),xi*(xi-1)*(2*eta+1);
                -4*xi*eta*(eta-1),   -2*(xi+1)*(xi-1)*(2*eta-1);
         -2*(2*xi+1)*(eta+1)*(eta-1),-4*xi*eta*(xi+1);
                -4*xi*eta*(eta+1),   -2*(xi+1)*(xi-1)*(2*eta+1);
         -2*(2*xi-1)*(eta+1)*(eta-1),-4*xi*eta*(xi-1);
                 8*xi*(eta^2-1),      8*eta*(xi^2-1)];
    end
    
    case 'Q16'
    % 16-Node Bicubic Quadrilateral
    if size(coord,2) < 2
      disp('Error two coordinates needed for the Q16 element')
    else
      r=coord(1); s=coord(2);     
shape(1) = -9.0/16.0*(r+1.0/3.0)*(r-1.0/3.0)*(r-1.0)*...
            (-9.0/16.0)*(s+1.0/3.0)*(s-1.0/3.0)*(s-1.0);
dhdr(1) = -9.0/16.0*((r+1.0/3.0)*(r-1.0/3.0)+(r+1.0/3.0)*(r-1.0)+...
         (r-1.0/3.0)*(r-1.0))*(-9.0/16.0)*(s+1.0/3.0)*(s-1.0/3.0)*(s-1.0);
dhds(1)= -9.0/16.0*(r+1.0/3.0)*(r-1.0/3.0)*(r-1.0)*(-9.0/16.0)*...
         ((s+1.0/3.0)*(s-1.0/3.0)+(s+1.0/3.0)*(s-1.0)+(s-1.0/3.0)*(s-1.0));
shape(5) = 27.0/16.0*(r+1.0)*(r-1.0/3.0)*(r-1.0)*...
           (-9.0/16.0)*(s+1.0/3.0)*(s-1.0/3.0)*(s-1.0);
dhdr(5) = 27.0/16.0*((r+1.0)*(r-1.0/3.0)+(r+1.0)*(r-1.0)+...
          (r-1.0/3.0)*(r-1.0))*(-9.0/16.0)*(s+1.0/3.0)*(s-1.0/3.0)*(s-1.0);
dhds(5) = 27.0/16.0*(r+1.0)*(r-1.0/3.0)*(r-1.0)*(-9.0/16.0)*...
         ((s+1.0/3.0)*(s-1.0/3.0)+(s+1.0/3.0)*(s-1.0)+(s-1.0/3.0)*(s-1.0));  
     
shape(6) = (-27.0/16.0)*(r+1.0)*(r+1.0/3.0)*(r-1.0)*...
            (-9.0/16.0)*(s+1.0/3.0)*(s-1.0/3.0)*(s-1.0);
dhdr(6) = (-27.0/16.0)*((r+1.0)*(r+1.0/3.0)+(r+1.0)*(r-1.0)+...
          (r+1.0/3.0)*(r-1.0))*(-9.0/16.0)*(s+1.0/3.0)*(s-1.0/3.0)*(s-1.0);
dhds(6) = (-27.0/16.0)*(r+1.0)*(r+1.0/3.0)*(r-1.0)*(-9.0/16.0)*...                    
        ((s+1.0/3.0)*(s-1.0/3.0)+ (s+1.0/3.0)*(s-1.0)+(s-1.0/3.0)*(s-1.0));
 
shape(2)= 9.0/16.0*(r+1.0)*(r+1.0/3.0)*(r-1.0/3.0)*...
            (-9.0/16.0)*(s+1.0/3.0)*(s-1.0/3.0)*(s-1.0);
dhdr(2) = 9.0/16.0*((r+1.0)*(r+1.0/3.0)+(r+1.0)*(r-1.0/3.0)+...
     (r+1.0/3.0)*(r-1.0/3.0))*(-9.0/16.0)*(s+1.0/3.0)*(s-1.0/3.0)*(s-1.0);
dhds(2) = 9.0/16.0*(r+1.0)*(r+1.0/3.0)*(r-1.0/3.0)*(-9.0/16.0)*...
         ((s+1.0/3.0)*(s-1.0/3.0)+(s+1.0/3.0)*(s-1.0)+(s-1.0/3.0)*(s-1.0));
 
shape(7) = 9.0/16.0*(r+1.0)*(r+1.0/3.0)*(r-1.0/3.0)*...
            (27.0/16.0)*(s+1.0)*(s-1.0/3.0)*(s-1.0);           
dhdr(7) = 9.0/16.0*((r+1.0)*(r+1.0/3.0)+(r+1.0)*(r-1.0/3.0)+...
         (r+1.0/3.0)*(r-1.0/3.0))*(27.0/16.0)*(s+1.0)*(s-1.0/3.0)*(s-1.0);
dhds(7) = 9.0/16.0*(r+1.0)*(r+1.0/3.0)*(r-1.0/3.0)*(27.0/16.0)*((s+1.0)*...
         (s-1.0/3.0)+(s+1.0)*(s-1.0)+(s-1.0/3.0)*(s-1.0));
shape(8) = 9.0/16.0*(r+1.0)*(r+1.0/3.0)*(r-1.0/3.0)*...
           (-27.0/16.0)*(s+1.0)*(s+1.0/3.0)*(s-1.0); 
dhdr(8) = 9.0/16.0*((r+1.0)*(r+1.0/3.0)+(r+1.0)*(r-1.0/3.0)+...
         (r+1.0/3.0)*(r-1.0/3.0))*(-27.0/16.0)*(s+1.0)*(s+1.0/3.0)*(s-1.0);
dhds(8) = 9.0/16.0*(r+1.0)*(r+1.0/3.0)*(r-1.0/3.0)*(-27.0/16.0)*...
             ((s+1.0)*(s+1.0/3.0)+(s+1.0)*(s-1.0)+(s+1.0/3.0)*(s-1.0));
shape(3) = 9.0/16.0*(r+1.0)*(r+1.0/3.0)*(r-1.0/3.0)*...
            9.0/16.0*(s+1.0)*(s+1.0/3.0)*(s-1.0/3.0);
dhdr(3) = 9.0/16.0*((r+1.0)*(r+1.0/3.0)+(r+1.0)*(r-1.0/3.0)+...
         (r+1.0/3.0)*(r-1.0/3.0))*9.0/16.0*(s+1.0)*(s+1.0/3.0)*(s-1.0/3.0);
dhds(3) = 9.0/16.0*(r+1.0)*(r+1.0/3.0)*(r-1.0/3.0)*9.0/16.0*...
         ((s+1.0)*(s+1.0/3.0)+(s+1.0)*(s-1.0/3.0)+(s+1.0/3.0)*(s-1.0/3.0));
shape(9) = (-27.0/16.0)*(r+1.0)*(r+1.0/3.0)*(r-1.0)*...
            9.0/16.0*(s+1.0)*(s+1.0/3.0)*(s-1.0/3.0);
dhdr(9) = (-27.0/16.0)*((r+1.0)*(r+1.0/3.0)+(r+1.0)*(r-1.0)+...
             (r+1.0/3.0)*(r-1.0))*9.0/16.0*(s+1.0)*(s+1.0/3.0)*(s-1.0/3.0);
dhds(9) = (-27.0/16.0)*(r+1.0)*(r+1.0/3.0)*(r-1.0)*9.0/16.0* ((s+1.0)*...
             (s+1.0/3.0)+(s+1.0)*(s-1.0/3.0)+(s+1.0/3.0)*(s-1.0/3.0));
         
shape(10) = 27.0/16.0*(r+1.0)*(r-1.0/3.0)*(r-1.0)*...
            9.0/16.0*(s+1.0)*(s+1.0/3.0)*(s-1.0/3.0);
dhdr(10) = 27.0/16.0*((r+1.0)*(r-1.0/3.0)+(r+1.0)*(r-1.0)+...
             (r-1.0/3.0)*(r-1.0))*9.0/16.0*(s+1.0)*(s+1.0/3.0)*(s-1.0/3.0);
dhds(10) = 27.0/16.0*(r+1.0)*(r-1.0/3.0)*(r-1.0)*9.0/16.0*((s+1.0)*...
                  (s+1.0/3.0)+(s+1.0)*(s-1.0/3.0)+(s+1.0/3.0)*(s-1.0/3.0));
shape(4) = -9.0/16.0*(r+1.0/3.0)*(r-1.0/3.0)*(r-1.0)*...
           9.0/16.0*(s+1.0)*(s+1.0/3.0)*(s-1.0/3.0);
dhdr(4) = -9.0/16.0*((r+1.0/3.0)*(r-1.0/3.0)+(r+1.0/3.0)*(r-1.0)+...
             (r-1.0/3.0)*(r-1.0))*9.0/16.0*(s+1.0)*(s+1.0/3.0)*(s-1.0/3.0);
dhds(4)= -9.0/16.0*(r+1.0/3.0)*(r-1.0/3.0)*(r-1.0)*9.0/16.0*...
         ((s+1.0)*(s+1.0/3.0)+(s+1.0)*(s-1.0/3.0)+(s+1.0/3.0)*(s-1.0/3.0));
 
shape(11) = -9.0/16.0*(r+1.0/3.0)*(r-1.0/3.0)*(r-1.0)*...
            (-27.0/16.0)*(s+1.0)*(s+1.0/3.0)*(s-1.0); 
dhdr(11) = -9.0/16.0*((r+1.0/3.0)*(r-1.0/3.0)+(r+1.0/3.0)*(r-1.0)+...
             (r-1.0/3.0)*(r-1.0))*(-27.0/16.0)*(s+1.0)*(s+1.0/3.0)*(s-1.0);
dhds(11) = -9.0/16.0*(r+1.0/3.0)*(r-1.0/3.0)*(r-1.0)*(-27.0/16.0)*...
                 ((s+1.0)*(s+1.0/3.0)+(s+1.0)*(s-1.0)+(s+1.0/3.0)*(s-1.0));
shape(12) = -9.0/16.0*(r+1.0/3.0)*(r-1.0/3.0)*(r-1.0)*...
           (27.0/16.0)*(s+1.0)*(s-1.0/3.0)*(s-1.0);
dhdr(12) = -9.0/16.0*((r+1.0/3.0)*(r-1.0/3.0)+(r+1.0/3.0)*(r-1.0)+...
              (r-1.0/3.0)*(r-1.0))*(27.0/16.0)*(s+1.0)*(s-1.0/3.0)*(s-1.0);
dhds(12) = -9.0/16.0*(r+1.0/3.0)*(r-1.0/3.0)*(r-1.0)*(27.0/16.0)*...
                 ((s+1.0)*(s-1.0/3.0)+(s+1.0)*(s-1.0)+(s-1.0/3.0)*(s-1.0));
shape(13) = 27.0/16.0*(r+1.0)*(r-1.0/3.0)*(r-1.0)*...
            (27.0/16.0)*(s+1.0)*(s-1.0/3.0)*(s-1.0); 
dhdr(13) = 27.0/16.0*((r+1.0)*(r-1.0/3.0)+(r+1.0)*(r-1.0)+...
              (r-1.0/3.0)*(r-1.0))*(27.0/16.0)*(s+1.0)*(s-1.0/3.0)*(s-1.0);
dhds(13) = 27.0/16.0*(r+1.0)*(r-1.0/3.0)*(r-1.0)*(27.0/16.0)*...
                 ((s+1.0)*(s-1.0/3.0)+(s+1.0)*(s-1.0)+(s-1.0/3.0)*(s-1.0));
shape(14) = (-27.0/16.0)*(r+1.0)*(r+1.0/3.0)*(r-1.0)*...
           (27.0/16.0)*(s+1.0)*(s-1.0/3.0)*(s-1.0);
dhdr(14) = (-27.0/16.0)*((r+1.0)*(r+1.0/3.0)+(r+1.0)*(r-1.0)+...
              (r+1.0/3.0)*(r-1.0))*(27.0/16.0)*(s+1.0)*(s-1.0/3.0)*(s-1.0);
dhds(14) = (-27.0/16.0)*(r+1.0)*(r+1.0/3.0)*(r-1.0)* (27.0/16.0)*...
                ((s+1.0)*(s-1.0/3.0)+ (s+1.0)*(s-1.0)+(s-1.0/3.0)*(s-1.0));
 
shape(15) = (-27.0/16.0)*(r+1.0)*(r+1.0/3.0)*(r-1.0)*...
            (-27.0/16.0)*(s+1.0)*(s+1.0/3.0)*(s-1.0);
dhdr(15) = (-27.0/16.0)*((r+1.0)*(r+1.0/3.0)+(r+1.0)*(r-1.0)+...
             (r+1.0/3.0)*(r-1.0))*(-27.0/16.0)*(s+1.0)*(s+1.0/3.0)*(s-1.0);
dhds(15) =  (-27.0/16.0)*(r+1.0)*(r+1.0/3.0)*(r-1.0)*(-27.0/16.0)*...
                 ((s+1.0)*(s+1.0/3.0)+(s+1.0)*(s-1.0)+(s+1.0/3.0)*(s-1.0));
shape(16) = 27.0/16.0*(r+1.0)*(r-1.0/3.0)*(r-1.0)*...
          (-27.0/16.0)*(s+1.0)*(s+1.0/3.0)*(s-1.0);
dhdr(16) = 27.0/16.0*((r+1.0)*(r-1.0/3.0)+(r+1.0)*(r-1.0)+...
             (r-1.0/3.0)*(r-1.0))*(-27.0/16.0)*(s+1.0)*(s+1.0/3.0)*(s-1.0);
dhds(16) = 27.0/16.0*(r+1.0)*(r-1.0/3.0)*(r-1.0)*(-27.0/16.0)*...
                 ((s+1.0)*(s+1.0/3.0)+(s+1.0)*(s-1.0)+(s+1.0/3.0)*(s-1.0));
N=shape(:);
dNdxi=[ dhdr(:) dhds(:)];
    end
    
   %% ==========================================================
   % 3D TETRAHEDRAL ELEMENTS
   % ==========================================================
   case 'H4'
    % 4-Node Linear Tetrahedron
    if size(coord,2) < 3
      disp('Error three coordinates needed for the H4 element')
    else
      xi=coord(1); eta=coord(2); zeta=coord(3);
      N=[1-xi-eta-zeta;
                    xi;
                   eta;
                  zeta];
      dNdxi=[-1  -1  -1;
              1   0   0;
              0   1   0;
              0   0   1];
    end
    
   case 'H10'
    % 10-Node Quadratic Tetrahedron (Added by Pranjal Saxena)
    if size(coord,2) < 3
        disp('Error: Three coordinates needed for the H10 element')
    else
        xi=coord(1); eta=coord(2); zeta=coord(3);
        N = [(1-xi-eta-zeta)*(1-2*xi-2*eta-2*zeta);
            xi*(2*xi-1);
            eta*(2*eta-1);
            zeta*(2*zeta-1);
            4*xi*(1-xi-eta-zeta);
            4*xi*eta;
            4*eta*(1-xi-eta-zeta);
            4*zeta*(1-xi-eta-zeta);
            4*xi*zeta;
            4*eta*zeta];
        dNdxi=[4*xi+4*eta+4*zeta-3,4*xi+4*eta+4*zeta-3,4*xi+4*eta+4*zeta-3;
            4*xi-1, 0, 0;
            0, 4*eta-1, 0;
            0, 0, 4*zeta-1;
            4-8*xi-4*eta-4*zeta, -4*xi, -4*xi;
            4*eta, 4*xi, 0;
            -4*eta,4-4*xi-8*eta-4*zeta, -4*eta;
            -4*zeta, -4*zeta,4-4*xi-4*eta-8*zeta;
            4*zeta, 0 , 4*xi;
            0, 4*zeta, 4*eta];
    end
    
    case 'H20'
    % 20-Node Cubic Tetrahedron (Added by Pranjal Saxena)
    if size(coord,2) < 3
        disp('Error: Three coordinates needed for the H20 element')
    else
        xi=coord(1); eta=coord(2); zeta=coord(3);
        N = [9*(1-xi-eta-zeta)*(1/3-xi-eta-zeta)*(2/3-xi-eta-zeta)/2;
            9*xi*(1/3-xi)*(2/3-xi)/2;
            9*eta*(1/3-eta)*(2/3-eta)/2;
            9*zeta*(1/3-zeta)*(2/3-zeta)/2;
            27*(1-xi-eta-zeta)*xi*(2/3-xi-eta-zeta)/2;
            27*(xi-1/3)*xi*(1-xi-eta-zeta)/2;
            27*xi*eta*(xi-1/3)/2;
            27*xi*eta*(eta-1/3)/2;
            27*eta*(eta-1/3)*(1-xi-eta-zeta)/2;
            27*eta*(1-xi-eta-zeta)*(2/3-xi-eta-zeta)/2;
            27*zeta*(1-xi-eta-zeta)*(2/3-xi-eta-zeta)/2;
            27*zeta*(zeta-1/3)*(1-xi-eta-zeta)/2;
            27*xi*zeta*(xi-1/3)/2;
            27*xi*zeta*(zeta-1/3)/2;
            27*eta*zeta*(zeta-1/3)/2;
            27*eta*zeta*(eta-1/3)/2;
            27*xi*eta*(1-xi-eta-zeta);
            27*zeta*eta*(1-xi-eta-zeta);
            27*xi*zeta*(1-xi-eta-zeta);
            27*xi*eta*zeta];
        dNdxi=[-(27*eta^2)/2 - 27*eta*xi - 27*eta*zeta + 18*eta - (27*xi^2)/2 - 27*xi*zeta + 18*xi - (27*zeta^2)/2 + 18*zeta-11/2, -(27*eta^2)/2 - 27*eta*xi - 27*eta*zeta + 18*eta - (27*xi^2)/2 - 27*xi*zeta + 18*xi - (27*zeta^2)/2 + 18*zeta - 11/2, - (27*eta^2)/2 - 27*eta*xi - 27*eta*zeta + 18*eta - (27*xi^2)/2 - 27*xi*zeta + 18*xi - (27*zeta^2)/2 + 18*zeta - 11/2;
            (27*xi^2)/2 - 9*xi + 1, 0, 0; 
            0, (27*eta^2)/2 - 9*eta + 1, 0;
            0, 0, (27*zeta^2)/2 - 9*zeta + 1;
            (27*eta^2)/2 + 54*eta*xi + 27*eta*zeta - (45*eta)/2 + (81*xi^2)/2 + 54*xi*zeta - 45*xi + (27*zeta^2)/2 - (45*zeta)/2 + 9, (9*xi*(6*eta + 6*xi + 6*zeta - 5))/2, (9*xi*(6*eta + 6*xi + 6*zeta - 5))/2;
            (9*eta)/2 + 36*xi + (9*zeta)/2 - 27*eta*xi - 27*xi*zeta - (81*xi^2)/2 - 9/2, -(xi*(27*xi - 9))/2, -(xi*(27*xi - 9))/2;
            (9*eta*(6*xi - 1))/2, (27*xi*(xi - 1/3))/2, 0;
            (27*eta*(eta - 1/3))/2, (9*xi*(6*eta - 1))/2, 0;
            -(27*eta*(eta - 1/3))/2, 36*eta + (9*xi)/2 + (9*zeta)/2 - 27*eta*xi - 27*eta*zeta - (81*eta^2)/2 - 9/2, -(27*eta*(eta - 1/3))/2;
            (9*eta*(6*eta + 6*xi + 6*zeta - 5))/2, (27*eta*(eta + xi + zeta - 1))/2 + (27*eta*(eta + xi + zeta - 2/3))/2 + (27*(eta + xi + zeta - 1)*(eta + xi + zeta - 2/3))/2, (9*eta*(6*eta + 6*xi + 6*zeta - 5))/2;
            (9*zeta*(6*eta + 6*xi + 6*zeta - 5))/2, (9*zeta*(6*eta + 6*xi + 6*zeta - 5))/2, (27*zeta*(eta + xi + zeta - 1))/2 + (27*zeta*(eta + xi + zeta - 2/3))/2 + (27*(eta + xi + zeta - 1)*(eta + xi + zeta - 2/3))/2;
            -(27*zeta*(zeta - 1/3))/2, -(27*zeta*(zeta - 1/3))/2, (9*eta)/2 + (9*xi)/2 + 36*zeta - 27*eta*zeta - 27*xi*zeta - (81*zeta^2)/2 - 9/2;
            (9*zeta*(6*xi - 1))/2, 0, (27*xi*(xi - 1/3))/2;
            (27*zeta*(zeta - 1/3))/2, 0, (9*xi*(6*zeta - 1))/2;
            0, (27*zeta*(zeta - 1/3))/2, (9*eta*(6*zeta - 1))/2;
            0, (9*zeta*(6*eta - 1))/2, (27*eta*(eta - 1/3))/2;
            -27*eta*(eta + 2*xi + zeta - 1), -27*xi*(2*eta + xi + zeta - 1), -27*eta*xi;
            -27*eta*zeta, -27*zeta*(2*eta + xi + zeta - 1), -27*eta*(eta + xi + 2*zeta - 1); 
            -27*zeta*(eta + 2*xi + zeta - 1), -27*xi*zeta, -27*xi*(eta + xi + 2*zeta - 1);
            27*eta*zeta, 27*xi*zeta, 27*eta*xi];
    end
    
   %% ==========================================================
   % 3D HEXAHEDRAL ELEMENTS
   % ==========================================================
   case 'H8'
    % 8-Node Linear Hexahedron
    if size(coord,2) < 3
      disp('Error three coordinates needed for the H8 element')
    else
      xi=coord(1); eta=coord(2); zeta=coord(3);
      I1=1/2-coord/2;
      I2=1/2+coord/2;
      N=[   I1(1)*I1(2)*I1(3);
            I2(1)*I1(2)*I1(3);
            I2(1)*I2(2)*I1(3);
            I1(1)*I2(2)*I1(3);
            I1(1)*I1(2)*I2(3);
            I2(1)*I1(2)*I2(3);
            I2(1)*I2(2)*I2(3);
            I1(1)*I2(2)*I2(3)   ];
      dNdxi=[   -1+eta+zeta-eta*zeta   -1+xi+zeta-xi*zeta  -1+xi+eta-xi*eta;
                 1-eta-zeta+eta*zeta   -1-xi+zeta+xi*zeta  -1-xi+eta+xi*eta;
                 1+eta-zeta-eta*zeta    1+xi-zeta-xi*zeta  -1-xi-eta-xi*eta;
                -1-eta+zeta+eta*zeta    1-xi-zeta+xi*zeta  -1+xi-eta+xi*eta;      
                -1+eta-zeta+eta*zeta   -1+xi-zeta+xi*zeta   1-xi-eta+xi*eta;
                 1-eta+zeta-eta*zeta   -1-xi-zeta-xi*zeta   1+xi-eta-xi*eta;
                 1+eta+zeta+eta*zeta    1+xi+zeta+xi*zeta   1+xi+eta+xi*eta;
                -1-eta-zeta-eta*zeta    1-xi+zeta-xi*zeta   1-xi+eta-xi*eta  ]/8;
    end
   
   case 'B20'
    disp(['Element ',type,' not yet supported'])
    if size(coord,2) < 3
      disp('Error three coordinates needed for the B20 element')
    else
      xi=coord(1); eta=coord(2); zeta=coord(3);
      N=[xi,eta,zeta];
      dNdxi=zeros(27,3);
    end
    
   case 'B27'
    disp(['Element ',type,' not yet supported'])
    if size(coord,2) < 3
      disp('Error three coordinates needed for the B27 element')
    else
      xi=coord(1); eta=coord(2); zeta=coord(3);
      N=zeros(27,1);
      dNdxi=zeros(27,3);
    end
    
   otherwise
    disp(['Element ',type,' not yet supported'])
    N=[]; dNdxi=[];
  end
end