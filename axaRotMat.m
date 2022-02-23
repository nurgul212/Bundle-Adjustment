function [R,dR,dRR] = axaRotMat(x)
%AXAROTMAT generates 3D axis-and-angle rotation matrix and its Jacobian.
%
%    R=axaRotMat(x) constructs the 3-by-3 Axis-Angle rotation matrix R by
%    a rotation angle and rotation axes.
%    R=cos(alpha)*I+(1-cos(alpha))*D_r+sin(alpha)*S_r, Equation(2.146)[1]
%     
%    [R,dR] = axaRotMat(x) also returns the Jacobian of R with respect to
%    the elements in input variable x. x is a 4-by-1 vector consists of a rotation angle 
%    x(1), and the 3-by-1 vector part x(2:4), i.e rotation axes.
%
%
%[1] Förstner, Wrobel, 2004, "Mathematical concepts in Photogrammetry", Ch. 2
%    of "Manual of Photogrammetry", McGlone et al., IAPRS.

if ~isvector(x) || any(size(x)~=[4,1])
    error('AXAROTMAT: Input variable must be 4-by-1');
end


alpha=x(1);        % a rotation angle
r=x(2:4);     % a rotation axes

I=eye(3,3);
D_r=r*r'; %symmetric dyadic matrix.
S_r= [0 -r(3) r(2);r(3) 0 -r(1); -r(2) r(1) 0];  % a skew symmetric matrix.
R=cos(alpha)*I+(1-cos(alpha))*D_r+sin(alpha)*S_r; 

if nargout>1
    % find Jacobian of R with respect to the angle alpha and axes r
    % respectively.
    s=sin(alpha);
    c=cos(alpha);
    
    dRda= [-s+r(1)^2*s, r(1)*r(2)*s-r(3)*c,r(1)*r(3)*s+r(2)*c;...
           r(1)*r(2)*s+r(3)*c,-s+r(2)^2*s,r(2)*r(3)*s-r(1)*c;...
           r(1)*r(3)*s-r(2)*c,r(2)*r(3)*s+r(1)*c,-s+r(3)^2*s]; 
   
    dRdr1=[2*r(1)*(1-c),r(2)*(1-c),r(3)*(1-c); r(2)*(1-c),0,-s;r(3)*(1-c),s,0]; 
    dRdr2=[0,r(1)*(1-c),s; r(1)*(1-c),2*r(2)*(1-c),r(3)*(1-c);-s,r(3)*(1-c),0]; 
    dRdr3=[0,-s,r(1)*(1-c);s,0,r(2)*(1-c);r(1)*(1-c),r(2)*(1-c),2*r(3)*(1-c)];  

    dR=[dRda(:),dRdr1(:),dRdr2(:),dRdr3(:)];
end

if nargout>2
    % Create a function that unroll the rotation matrix.
    f=@(x)reshape(axaRotMat(x),9,1);
    dRR=jacapprox(f,x,1e-6,{});
end
