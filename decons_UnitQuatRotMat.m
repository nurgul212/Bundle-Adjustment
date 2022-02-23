function x = decons_UnitQuatRotMat(R )
%DECONS_UNITQUATRODMAT returns unit quaternion parameters as a 3-by-1 vector
% from deconstructing a rotation matrix R.

% derive a scalar s and vector part v by deconstruction of R
s=(1/2)*sqrt(1+R(1,1)+R(2,2)+R(3,3));
v(1)=(1/(4*s))*(R(3,2)-R(2,3));
v(2)=(1/(4*s))*(R(1,3)-R(3,1));
v(3)=(1/(4*s))*(R(2,1)-R(1,2));
x=[s;v']; 
 
end


