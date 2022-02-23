function [R,dR,dRR] = UnitQuatRotMat(x)
%UNITQUATROTMAT generates 3D rotation matrix by unit quaternion and its Jacobian.
%
%    R = UnitQuatRotMat(x) constucts the 3-by-3 rotation matrix by
%    unit quaternion. 
%    R= ((s*s-v'*v)*I+2*D_v+2*s*S_v).
%    
%    [R,dR]= UnitQuatRotMat(x) also returns the Jacobian of R with respect to
%    the each elements of unit quaternion q=[s,v]', where s is a scalar part 
%    v is a vector part.
%   


if ~isvector(x) || any(size(x)~=[4,1])
    error('UnitQUATROTMAT: Input variable must be 4-by-1')
end


s=x(1);   % scalar part
v=x(2:4); % a vector part

I=eye(3,3);
D_v=v*v';     %symmetric dyadic matrix.
S_v= [0 -v(3) v(2);v(3) 0 -v(1); -v(2) v(1) 0];  % a skew symmetric matrix.
R=((s*s-v'*v)*I+2*D_v+2*s*S_v);  % Rotation by Unit quaternion 

if nargout>1
    % find Jacobian of R with respect to a scalar s and vector v respectively.
    dRds=2*[s,-v(3),v(2);v(3),s,-v(1);-v(2),v(1),s]; 
    dRdv1=2*[v(1),v(2),v(3);v(2),-v(1),-s;v(3),s,-v(1)];  
    dRdv2=2*[-v(2),v(1),s;v(1),v(2),v(3);-s,v(3),-v(2)];  
    dRdv3=2*[-v(3),-s,v(1);s,-v(3),v(2);v(1),v(2),v(3)];  
    
    dR=[dRds(:),dRdv1(:),dRdv2(:),dRdv3(:)];
end

if nargout>2
    % Create a function that unroll the rotation matrix.
    f=@(x)reshape(UnitQuatRotMat(x),9,1);
    dRR=jacapprox(f,x,1e-6,{});
end
