function plot_mesh(X,connect,elem_type,se)

% function plot_mesh(X,connect,elem_type,linespec)
% 
% plots a nodal mesh and an associated connectivity.  X is
% the nodal coordinates, connect is the connectivity, and
% elem_type is either 'L2', 'L3', 'T3', 'T6', 'Q4', or 'Q9' 
% depending on the element topology.
  
if ( nargin < 4 )
   se='wo-';
end

holdState=ishold;
hold on

% fill X if needed
if (size(X,2) < 3)
   for c=size(X,2)+1:3
      X(:,c)=[zeros(size(X,1),1)];
   end
end

for e=1:size(connect,1)
  
   if ( strcmp(elem_type,'Q9') )       % 9-node quad element
      ord=[1,2,3,6,9,8,7,4,1];
   elseif ( strcmp(elem_type,'Q8') )  % 8-node quad element
 %      ord=[1,5,2,6,3,7,4,8,1];
      ord=[1,2,3,4,1];
   elseif ( strcmp(elem_type,'T3') )  % 3-node triangle element
      ord=[1,2,3,1];
   elseif ( strcmp(elem_type,'T4') )  % 3-node triangle element
      ord=[1,2,3,1];   
   elseif ( strcmp(elem_type,'T6') )  % 6-node triangle element
      ord=[1,4,2,5,3,6,1];
   elseif ( strcmp(elem_type,'T10') )  % 10-node triangle element
      ord=[1,4,5,2,6,7,3,8,9,1];
   elseif ( strcmp(elem_type,'T15') )  % 15-node triangle element
      ord=[1,4,5,6,2,7,8,9,3,10,11,12,1];
   elseif ( strcmp(elem_type,'Q4') )  % 4-node quadrilateral element
      ord=[1,2,3,4,1];
   elseif ( strcmp(elem_type,'L2') )  % 2-node line element
      ord=[1,2];   
   elseif ( strcmp(elem_type,'L3') )  % 3-node line element
      ord=[1,3,2];   
   elseif ( strcmp(elem_type,'H4') )  % 4-node tet element
      ord=[1,2,4,1,3,4,2,3];   
   elseif ( strcmp(elem_type,'H4b') )  % 4-node tet element
      ord=[1,2,4,1,3,4,2,3];  
   elseif ( strcmp(elem_type,'H10') )  % ps %10-node tet element
      ord=[1,5,2,9,4,8,1,7,3,10,4,8,1,5,2,6,3];
   elseif ( strcmp(elem_type,'H20') )  % ps %10-node tet element
      ord=[1,5,6,2,13,14,4,12,11,1,10,9,3,16,15,4,12,11,1,5,6,2,7,8,3];
   elseif ( strcmp(elem_type,'H8') )  % 8-node hex element
      ord=[1,5,6,2,3,7,8,4,1,2,3,4,8,5,6,7];  
   elseif ( strcmp(elem_type,'B20') )  % ps %20-node brick element
      ord=[1,9,2,10,3,11,4,12,1,17,5,13,6,14,7,15,8,16,5,13,6,18,2,10,3,19,7,15,8,20,4];
   elseif ( strcmp(elem_type,'Q25') )  % 8-node brick element
      ord=[1,2 3 4 5 10 15 20 25 24 23 22 21 16 11 6 1];   
   end
   
   for n=1:size(ord,2)
      xpt(n)=X(connect(e,ord(n)),1);
      ypt(n)=X(connect(e,ord(n)),2);      
      zpt(n)=X(connect(e,ord(n)),3);
   end
   plot3(xpt,ypt,zpt,se,LineWidth=1.5)
end

rotate3d on
axis equal
      
if ( ~holdState )
  hold off
end
