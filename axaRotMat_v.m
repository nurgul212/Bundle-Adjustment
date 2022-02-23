function [v,J]=axaRotMat_v(varargin)
%AXAROTMAT_V Vector version of axaRotMat. Used for testing only.
%
%   v=axaRotMat_v(x) returns the output of axaRotMat as a column
%   vector.
%
%   [v,J]=axaRotMat_v('selftest',...) performs a selftest of the
%   axaRotMat function.

if ~isempty(varargin)&& ischar(varargin{1})&& strcmp(varargin{1},'selftest')
    
    [v,J]=selftest(varargin{2:end});
else
    % Call axaRotMat and unroll the result.
    [R,J]=axaRotMat(varargin{:});
    v=R(:);
end

function [Ja,Jn]=selftest(varargin)

        x=[1,3,5,7]'; 
  
    if ( norm(x(2:4))~=1)
    x=[x(1);x(2:4)/(norm(x(2:4)))];  % normalize the vector part
    end

%end Numerical Jacobian.
Jn=jacapprox(mfilename,x,1e-6,{});
% Analytical Jacobian. 
[~,Ja]=axaRotMat_v(x);
end
end