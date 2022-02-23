function [C,J,JJc] = Rdcm( x)
%RDCM generates the 3-by-3 reduced direction cosine matrix and its Jacobian
  

if ~isvector(x) || any(size(x)~=[6,1])
    error('RDCM: Input variable must be 6-by-1');
end

r1=x(1:3);
r2=x(4:6);
r3=cross(r1,r2);
C=[r1,r2,r3];

J_c11=[1 0 0;0 0 -r2(3);0 0 r2(2)];
J_c12=[0 0 r2(3);1 0 0;0 0 -r2(1)];
J_c13=[0 0 -r2(2);0 0 r2(1);1 0 0];

J_c21=[0 1 0;0 0 r1(3);0 0 -r1(2)];
J_c22=[0 0 -r1(3);0 1 0;0 0 r1(1)];
J_c23=[0 0 r1(2);0 0 -r1(1);0 1 0];

J=[J_c11(:),J_c12(:),J_c13(:),J_c21(:),J_c22(:),J_c23(:)];

if nargout>2
    % Numerical approximation.
    f=@(x)reshape(Rdcm(x),9,1);
    JJc=jacapprox(f,x,1e-6,{});
end


