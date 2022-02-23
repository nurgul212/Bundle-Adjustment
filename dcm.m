function [C,J, JJc]=dcm(x)
%RDCM generates 3-by-3 direction cosine matrix and its Jacobian.
  

if ~isvector(x) || any(size(x)~=[9,1])
    error('DCM: Input variable must be 9-by-1');
end

r1=x(1:3);
r2=x(4:6);
r3=x(7:9);
C=[r1,r2,r3];

J=eye(9,9);

if nargout>2
    % Numerical approximation.
    f=@(x)reshape(dcm(x),9,1);
    JJc=jacapprox(f,x,1e-6,{});
end
