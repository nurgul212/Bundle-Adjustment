function [c,Jc,JJc]=dcm_c(x)
%DCM_C Constraint function for direction cosine matrix (DCM) representation of a rotation.
%
%   c=DCM_C(X,...) returns the constraints that need to be satisfied for
%   the DCM representation of a rotation. i.e. Ci'Ci-1=0; CiCj=0 (j~=i);
%   
%
%   [c,Jc]=... also returns the 6-by-9 Jacobian of c w.r.t. each element.
%   Jacobian=[2C1'  0     0  
%             0   2C2'   0
%             0    0    2C3'
%            C2'  C1'    0
%            C3'   0    C1'
%            0    C3'   C2' ]


if ~isvector(x) || any(size(x)~=[9,1])
    error('DCM_C: x must be 9-by-1 vector');
end

C=[x(1:3),x(4:6),x(7:end)]; %DCM
c=zeros(6,1);

for i=1:3
    c(i)=C(:,i)'*C(:,i)-1;
end  
c(4)=C(:,1)'*C(:,2);
c(5)=C(:,1)'*C(:,3);
c(6)=C(:,2)'*C(:,3);

if nargout>1
    % Find Jacobian of constraints c
    Jc1=[2*C(:,1)';zeros(1,3);zeros(1,3);C(:,2)';C(:,3)';zeros(1,3)];
    Jc2=[zeros(1,3);2*C(:,2)';zeros(1,3);C(:,1)';zeros(1,3);C(:,3)'];
    Jc3=[zeros(1,3);zeros(1,3);2*C(:,3)';zeros(1,3);C(:,1)';C(:,2)'];
    Jc=[Jc1,Jc2,Jc3];
end
 
if nargout>2
    % Numerical approximation.
    JJc=jacapprox(mfilename,x,1e-6,{});
end
    

