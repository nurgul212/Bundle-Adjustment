function [ c,Jc,JJc ] = UnitQuat_c(x)
%UNITQUAT_C Constraint function for unit quaternion representation of a rotation.
%
%   c=UNITQUAT_C(X,...) returns the constraints that need to be satisfied for
%   the unit quaternion representation of a rotation. X=[s;v], where s is
%   the scalar part and v is the 3-by-1 vector 
%
%   [c,Jc]=... also returns the Jacobian of c with respect to x.

if ~isvector(x) || any(size(x)~=[4,1])
    error('UNITQUAT_C: x must be 4-by-1 vector');
end
c=x'*x-1;  % s^2+v(1)^2+v(2)^2+v(3)^2
if nargout>1
% s=x(1);   %  scalar part
% v=x(2:4); %  vector part
Jc=2*x';
end

if nargout>2
    % Numerical approximation.
    JJc=jacapprox(mfilename,x,1e-6,{});
end
end

