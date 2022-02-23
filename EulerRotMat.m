function [R,dR] = EulerRotMat(angles,seq,isFixed)
%EULERROTMAT Generate 3D Euler angle rotation matrix and its Jacobian.
%
%   R1=EULERROTMAT(ANG,SEQ) constructs the 3-by-3 Euler rotation matrix R1
%   by the angles in the 3-vector ANG=[a,b,c] and with the angle sequence
%   specified in 3-vector SEQ. The returned rotation matrix R1 correponds to
%   rotation about fixed reference axes [1]. SEQ can be any 3-letter
%   sequence of the letters 'x', 'y', and 'z'.
%
%   R2=EULERROTMAT(ANG,SEQ,TRUE) constructs the 3-by-3 Euler rotation matrix
%   that corresponds to rotation about moving axes [1].
%
%   [R,dR]=... also returns the Jacobian of R with respect to the elements
%   in ANG. The ordering of the rows of dR correspond to a columnwise
%   unrolling of R, i.e. R(:). The partial derivatives are computed by the
%   method described in [2].
%
%[1] Förstner, Wrobel, 2004, "Mathematical concepts in Photogrammetry", Ch. 2
%    of "Manual of Photogrammetry", McGlone et al., IAPRS.
%[2] Lucas, 1963, "Differentiation of the Orientation Matrix by Matrix
%    Multipliers", Photogrammetric Engineering 29(4):708-715.

if nargin<3, isFixed=true; end

if ~isvector(angles) || any(size(angles)~=[3,1])
    error('EULERROTMAT: Angles must be 3-by-1');
end
if ~isvector(seq) || ...
        any(size(seq)~=[1,3])|| ...
        ~ischar(seq) || ...
        ~(strcmp(seq,'xyz') || strcmp(seq,'ats') || strcmp(seq,'zxz'))
 error('EULERROTMAT: Axis sequence must be char 1-by-3-vector with chars (x, y, z) or (a,t,s)');
end

switch seq
  case 'xyz' 
    if isFixed
        % Rotation about fixed axes, so the multiplication sequence is
        % M=Mz(c)My(b)Mx(a). The sign change on the angles and the Jacobian
        % is due to Lucas' definition of positive angles as clockwise.
        if nargout<2
            R=lucas(-flipud(angles),fliplr(seq));
        else
            [R,dRda,dRdb,dRdc]=lucas(-flipud(angles),fliplr(seq));
            dR=-[dRdc(:),dRdb(:),dRda(:)];
        end
    else
        % Rotation about moving axes, so the multiplication sequence is
        % M=Mx(a)My(b)Mz(c). The sign change on the angles and the
        % Jacobian is due to Lucas' definition of positive angles as
        % clockwise.
        if nargout<2
            R=lucas(-angles,seq);
        else
            [R,dRda,dRdb,dRdc]=lucas(-angles,seq);
            dR=-[dRda(:),dRdb(:),dRdc(:)];
        end
    end
    
 case 'zxz'    
      if isFixed
        % Rotation about fixed axes, so the multiplication sequence is
        % M=Mz(c)Mx(b)Mz(a). The sign change on the angles and the Jacobian
        % is due to Lucas' definition of positive angles as clockwise.
        if nargout<2
            R=lucas(-flipud(angles),seq);
        else
            [R,dRda,dRdb,dRdc]=lucas(-flipud(angles),seq);
            dR=-[dRdc(:),dRdb(:),dRda(:)];
        end   
      else  % rotation as transformation
         if nargout<2
            R=lucas(-angles,seq);
        else
            [R,dRda,dRdb,dRdc]=lucas(-angles,seq);
            dR=-[dRdc(:),dRdb(:),dRda(:)];
         end 
      end
      
case 'ats'  
    Angles=angles-[0,0,pi]';  % azimuth, tilt, swing angles are clockwise
    angles=[-Angles(1);Angles(2);Angles(3)];
     seq='zxz';
      if isFixed
        % Rotation about fixed axes, so the multiplication sequence is
        % M=Mz(c)Mx(b)Mz(a). The sign change on the angles and the Jacobian
        % is due to Lucas' definition of positive angles as clockwise.
        if nargout<2
            R=lucas(flipud(angles),seq);
        else
            [R,dRda,dRdb,dRdc]=lucas(-flipud(angles),seq);
            dR=-[dRdc(:),dRdb(:),dRda(:)];
        end   
      else  % rotation as transformation
         if nargout<2
            R=lucas(angles,seq);
        else
            [R,dRda,dRdb,dRdc]=lucas(angles,seq);
            dR=[dRdc(:),dRdb(:),dRda(:)];
         end 
      end 
      
    otherwise
    error('EulerRotMat:Illegal axis sequence')
end
end       
