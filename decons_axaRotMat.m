function x = decons_axaRotMat( R )
%DECONS_AXAROTMAT returns axis-angle parameters from deconstructing a
%rotation matrix and errer between axis-angle used for constructing and
%deconstructing rotation

 % construct a rotation matrix R by the given parameters.   

 a=-[R(2,3)-R(3,2);R(3,1)-R(1,3);R(1,2)-R(2,1)];
 alpha = atan2(norm(a),trace(R)-1);   % derive a rotation angle from R
 r=a/norm(a);                         % derive rotation axes from R
    
 x=[alpha;r];     % find rotation angle and axes from deconstructing R 
     
end