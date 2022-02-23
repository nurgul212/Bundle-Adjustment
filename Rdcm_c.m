function [c,Jc,JJc] = Rdcm_c(x)
%RDCM_C Constraint function for reduced direction cosine matrix (RDCM) 
% representation of a rotation.
%
%
%   c=RDCM_C(X,...) returns the constraints that need to be satisfied for
%   the RDCM representation of a rotation. i.e. Ci'Ci-1=0; CiCj=0 (j~=i);
%   
%
%   [c,Jc]=... also returns the 3-by-6 Jacobian of c w.r.t. C1,C2.
%   Jacobian=[ 2C1'  0       
%               0   2C2'   
%               C2'  C1' ]   
%           

if ~isvector(x) || any(size(x)~=[6,1])
    error('RDCM_C: Input variable must be 6-by-1');
end

C=[x(1:3),x(4:end),cross(x(1:3),x(4:end))];  %RDCM    
c=zeros(3,1);
for i=1:2
    c(i)=C(:,i)'*C(:,i)-1;
end
c(3)=C(:,1)'*C(:,2);

if nargout>1
    % Find Jacobian of constraints c
    Jc1=[2*C(:,1)';zeros(1,3);C(:,2)'];
    Jc2=[zeros(1,3);2*C(:,2)';C(:,1)'];
    Jc=[Jc1,Jc2];
end

if nargout>2
    % Numerical approximation.
    JJc=jacapprox(mfilename,x,1e-6,{});
end
    


