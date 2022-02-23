function [v,J]=lucas_v(varargin)
%LUCAS_V Vector version of Lucas. Used for testing only.
%
%   [v,J]=lucas_v(angles,seq) returns the output of Lucas as a column
%   vector.
%
%   [v,J]=lucas_v('selftest',...) performs a selftest of the lucas function.

if ~isempty(varargin) && ischar(varargin{1}) && strcmp(varargin{1},'selftest')
    [v,J]=selftest(varargin{2:end});
else
    [angles,seq]=deal(varargin{:});
    % Call Lucas and unroll the result.
    [M,dMda,dMdb,dMdc]=lucas(angles,seq);
    v=M(:);
    J=[dMda(:),dMdb(:),dMdc(:)];
end

function [Ja,Jn]=selftest(angles,seq)

if nargin<1, angles=[3,5,7]'*pi/180; end
if nargin<2, seq='xyz'; end

% Numerical Jacobian.
Jn=jacapprox(mfilename,angles,1e-6,{seq});

% Analytical Jacobian. 
[~,Ja]=lucas_v(angles,seq);
