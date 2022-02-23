function [v,J]=EulerRotMat_v(varargin)
%EULERROTMAT_V Vector version of EulerRotMat. Used for testing only.
%
%   v=EulerRotMat_v(ANG,SEQ) returns the output of EulerRotMat as a column
%   vector.
%
%   [v,J]=EulerRotMat_v('selftest',...) performs a selftest of the
%   EulerRotMat function.

if ~isempty(varargin) && ischar(varargin{1}) && strcmp(varargin{1},'selftest')
    [v,J]=selftest(varargin{2:end});
else
    % Call Lucas and unroll the result.
    [R,J]=EulerRotMat(varargin{:});
    v=R(:);
end

function [Ja,Jn]=selftest(varargin)

if length(varargin)<1
    angles=[3,5,7]'*pi/180;
else
    angles=varargin{1};
end

if length(varargin)<2
    seq='xyz';
else
    seq=varargin{2};
end

if length(varargin)<3
    isFixed=true;
else
    isFixed=varargin{3};
end

% Numerical Jacobian.
Jn=jacapprox(mfilename,angles,1e-6,{seq,isFixed});

% Analytical Jacobian. 
[~,Ja]=EulerRotMat_v(angles,seq,isFixed);
end
end