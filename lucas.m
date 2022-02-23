function [M,dMda,dMdb,dMdc]=lucas(angles,seq)
%LUCAS Construct rotation matrix and Jacobian according to Lucas (1963).
%
%   M=lucas(ANG,SEQ), where ANG=[a,b,c] are rotation angles and SEQ=[i,j,k]
%   is an axis sequence, computes the 3-by-3 rotation matrix M that is
%   formed according to [1], i.e. as
%
%          M = M_i(a) M_j(b) M_k(c),
%
%   where i, j, k are 'x', 'y', or 'z', and M_{i,j,k}(v) are the planar
%   rotation matrices about the corresponding axes. A positive rotation
%   angle v corresponds to a clockwise rotation.
%
%   [M,dMda,dMdb,dMdc]=... also returns the 3-by-3 partial derivative
%   matrices according to [1].
%
%References:
%[1] Lucas 1963, "Differentiation of the Orientation Matrix by Matrix
%    Multipliers", Photogrammetric Engineering 29(4), 708-715.

% if ~isvector(angles) || any(size(angles)~=[3,1])
%     error('LUCAS: Angles must be 3-by-1');
% end
% 
if ~isvector(seq) || any(size(seq)~=[1,3]) || ~ischar(seq) || ...
        any(~ismember(lower(seq),'xyz'))
  
    error('LUCAS: Axis sequence must be char 1-by-3-vector with chars x, y, z');
end


% Matrices defined by eqn. (6).
A=planrot(angles(1),seq(1));
B=planrot(angles(2),seq(2));
C=planrot(angles(3),seq(3));
M=A*B*C;
    
% Eqn. (7).
dMda=Pi(seq(1))*M;
% Eqn. (8).
dMdb=A*B*Pi(seq(2))*C;
% Eqn. (9).
dMdc=M*Pi(seq(3));


function M=planrot(a,i)
% Planar rotations defined by eqn. (1) in [1]. The rotation angle a is
% given in radians, i='x',' y'', 'z' is the rotation axis.

M=eye(3);

% 2D rotation.
R=[cos(a),sin(a);-sin(a),cos(a)];

% Place it into the proper place in M.
switch i
  case 'x'
    % Rotation about the x axis.
    M(2:3,2:3)=R;
  case 'y'
    M([1,3],[1,3])=R';
  case 'z'
    M(1:2,1:2)=R;
  otherwise
    error('Bad rotation axis.')
end

function P=Pi(i)
% Skew-symmetric helper matrices defined by eqn. (3) of [1].

P=zeros(3);

% Non-zero block.
R=[0,1;-1,0];
switch i
  case 'x'
    P(2:3,2:3)=R;
  case 'y'
    P([1,3],[1,3])=R';
  case 'z'
    P(1:2,1:2)=R;
  otherwise
    error('Bad axis.')
end

