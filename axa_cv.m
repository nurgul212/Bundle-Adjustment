function [c,A]=axa_cv(varargin)
%AXA_CV Vector version of AXA_C. Used for testing only.
%
%   c=axa_cv(x,...) returns the output of AXA_CV as a column vector.
%
%   [c,A]=axa_cv('selftest',...) performs a selftest of the axa_c function.

if ~isempty(varargin) && ischar(varargin{1}) && strcmp(varargin{1},'selftest')
    [c,A]=selftest(varargin{2:end});
else
    % Call AXA_C and unroll the result.
    [c,A]=axa_c(varargin{:});
end

function [Ja,Jn]=selftest(x,varargin)

if nargin<1, x=[1,2,3,4]'; 
end

% Numerical Jacobian.
Jn=jacapprox(mfilename,x,1e-6);

% Analytical Jacobian. 
[~,Ja]=axa_cv(x);
