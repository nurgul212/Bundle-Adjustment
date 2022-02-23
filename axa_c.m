function [c,Jc,JJc]=axa_c(x)
%AXA_C Constraint function for axis-and-angle representation of a rotation.
%
%   c=AXA_C(X,...) returns the constraints that need to be satisfied for
%   the axis-and-angle representation of a rotation. X=[alpha;r], where
%   the scalar alpha is the rotation angle and the 3-by-1 vector r is the
%   unit vector the corresponds to the rotation axis.
%
%   [c,Jc]=... also returns the Jacobian of c w.r.t. X.

if ~isvector(x) || any(size(x)~=[4,1])
    error('AXA_C: x must be 4-by-1 vector');
end

%alpha=x(1);
r=x(2:4);

c=r'*r-1;  % r(1)^2+r(2)^2+r(3)^2

if nargout>1
    Jc=[0,2*r'];
end

if nargout>2
    % Numerical approximation.
    JJc=jacapprox(mfilename,x,1e-6,{});
end
