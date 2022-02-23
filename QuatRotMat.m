function [R,dR,dRR] = QuatRotMat( x )
%QUATROTMAT generates 3D rotation matrix by quaternion.
%
%    R = QuatRotMat( x ) constucts the 3-by-3 rotation matrix by
%    quaternion. 
%    R= (1/(s*s+v(1)^2+v(2)^2+v(3)^2))*((s*s-v'*v)*I+2*D_v+2*s*S_v),
%    Equation(2.169). [1]
%
%    
%    [R,dR]= QuatRotMat(x) also returns the Jacobian of R with respect to
%    the each elements of  quaternion q=[s,v]', where s is a scalar part 
%    v is a vector part.
%    
%
%[1] Förstner, Wrobel, 2004, "Mathematical concepts in Photogrammetry", Ch. 2
%    of "Manual of Photogrammetry", McGlone et al., IAPRS.

if ~isvector(x) || any(size(x)~=[4,1])
    error('QUATROTMAT: Input variable must be 4-by-1');
end

s=x(1);    % a scalar part
v1=x(2);    
v2=x(3);
v3=x(4);
v=[v1,v2,v3]';% 3-by-1 vector part 
    
I=eye(3,3);
D_v=v*v';    %symmetric dyadic matrix.
S_v= [0 -v3 v2;v3 0 -v1; -v2 v1 0];  % a skew symmetric matrix.
R=(1/(s*s+v1^2+v2^2+v3^2))*((s*s-v'*v)*I+2*D_v+2*s*S_v);   % Rotation by quaternion

if nargout>1
    % find Jacobian of R with respect to the quaternian parameters s, v1, v2,
    % and v3 respectively.

    n=(s^2+v1^2+v2^2+v3^2);
    a=s^2-v1^2-v2^2-v3^2;
    b=s^2-v1^2+v2^2+v3^2;
    c=s^2+v1^2-v2^2+v3^2;
    d=s^2+v1^2+v2^2-v3^2;
    
    dRds=1/n^2*([4*s*(v2^2+v3^2),   2*v3*a-4*s*v1*v2, -2*v2*a-4*s*v1*v3;...
                 -2*v3*a-4*s*v1*v2,  4*s*(v1^2+v3^2),  2*v1*a-4*s*v2*v3;...
                 2*v2*a-4*s*v1*v3, -2*v1*a-4*s*v2*v3, 4*s*(v1^2+v2^2)]);
    
    dRdv1=1/n^2*([4*v1*(v2^2+v3^2), 2*v2*b+4*s*v1*v3, 2*v3*b-4*s*v1*v2;...
                  2*v2*b-4*s*v1*v3, -4*v1*(s^2+v2^2), -2*s*b-4*v1*v2*v3;...
                  2*v3*b+4*s*v1*v2,  2*s*b-4*v1*v2*v3, -4*v1*(s^2+v3^2)]);
    
    dRdv2=1/n^2*([-4*v2*(s^2+v1^2), 2*v1*c+4*s*v2*v3,  2*s*c-4*v1*v2*v3;...
                  2*v1*c-4*s*v2*v3,  4*v2*(v1^2+v3^2),  2*v3*c+4*s*v1*v2;...
                  -2*s*c-4*v1*v2*v3, 2*v3*c-4*s*v1*v2,  -4*v2*(s^2+v3^2)]);
    
    dRdv3=1/n^2*([-4*v3*(s^2+v1^2), -2*s*d-4*v1*v2*v3, +2*v1*d-4*s*v2*v3;...
                  2*s*d-4*v1*v2*v3,  -4*v3*(s^2+v2^2),  2*v2*d+4*s*v1*v3;...
                  2*v1*d+4*s*v2*v3,  2*v2*d-4*s*v1*v3,  4*v3*(v1^2+v2^2)]);
       
    dR=[dRds(:),dRdv1(:),dRdv2(:),dRdv3(:)];
end

if nargout>2
    % Create a function that unroll the rotation matrix.
    f=@(x)reshape(QuatRotMat(x),9,1);
    dRR=jacapprox(f,x,1e-6,{});
end
