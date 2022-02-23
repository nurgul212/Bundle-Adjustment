function [R,dR,dRR] = RodriRotMat( m)
%RODRIROTMAT generates  a 3D rotation matrix by Rodriguez representation
%and its Jacobian.
%
%    R = RodriRotMat( m ) constucts the 3-by-3 rotation matrix by
%    Rodriguez representation, i.e. R=(1/(norm(m))^2)*((4-m'*m)*I+2*D_m+4*S_m);
%    Equation(2.174),[1].
%    q=[1,a/2,b/2,c/2]', where m=[a,b,c]'.
%
%    [R,dR]=RodriRotMat( m ) also returns the Jacobian of R with respect to a,b and c elements of m.
%    input variable m is a 3-by-1 vector.
%
%
%
%[1] Förstner, Wrobel, 2004, "Mathematical concepts in Photogrammetry", Ch. 2
%    of "Manual of Photogrammetry", McGlone et al., IAPRS.

if ~isvector(m) || any(size(m)~=[3,1])
    error('RODRIROTMAT: Input variable must be 3-by-1');
end

a=m(1);   
b=m(2); 
c=m(3);

I=eye(3,3);
D_m=[a;b;c]*[a;b;c]'; %symmetric dyadic matrix.
S_m= [0 -c b;c 0 -a; -b a 0];  % a skew symmetric matrix.
R=(1/(4+a*a+b*b+c*c))*((4-[a;b;c]'*[a;b;c])*I+2*D_m+4*S_m);  % rotation by Rodriguez


if nargout>1
    % find Jacobian of R with respect to a,b and c. 

    dRda= 1/(4+m'*m)^2.*[4*a*(b*b+c*c),2*b*(-a*a+b*b+c*c+4)+8*a*c, 2*c*(-a*a+b*b+c*c+4)-8*a*b;...
        2*b*(-a*a+b*b+c*c+4)-8*a*c,-4*a*(b*b+4),-4*(-a*a+b*b+c*c+a*b*c+4);...
        2*c*(-a*a+b*b+c*c+4)+8*a*b,4*(-a*a+b*b+c*c-a*b*c+4),-4*a*(c*c+4)];
    
    dRdb= 1/(4+m'*m)^2.*[-4*b*(a*a+4),2*a*(a*a-b*b+c*c+4)+8*b*c, 4*(a*a-b*b+c*c-a*b*c+4);...
        2*a*(a*a-b*b+c*c+4)-8*b*c,4*b*(a*a+c*c),2*c*(a*a-b*b+c*c+4)+8*a*b;...
        -4*(a*a-b*b+c*c+a*b*c+4),2*c*(a*a-b*b+c*c+4)-8*a*b,-4*b*(c*c+4)];

    
    dRdc= 1/(4+m'*m)^2.*[-4*c*(a*a+4),-4*(a*a+b*b-c*c+a*b*c+4), 2*a*(a*a+b*b-c*c+4)-8*b*c;...
        4*(a*a+b*b-c*c-a*b*c+4),-4*c*(b*b+4),2*b*(a*a+b*b-c*c+4)+8*a*c;...
        2*a*(a*a+b*b-c*c+4)+8*b*c,2*b*(a*a+b*b-c*c+4)-8*a*c,4*c*(a*a+b*b)];


   dR=[dRda(:),dRdb(:),dRdc(:)];
end

if nargout>2
    % Create a function that unroll the rotation matrix.
    f=@(x)reshape(RodriRotMat(x),9,1);
    dRR=jacapprox(f,m,1e-6,{});
end
