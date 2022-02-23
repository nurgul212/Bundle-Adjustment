function m= decons_RodriRotMat(R)
%DECONS_RODRIROTMAT deconstructs a rotation matrix R to return Rodriguez
%parameter m as a 3-by-1 vector
%
%
%   m= decons_RodriRotMat( R ) returns Rodriguez parameter m by using the 
%   relation between axis-angle representation and Rodriguez, i.e. rotation 
%   axes r=m/norm(m), so that m=r*norm(m). From e.g (2.175) [1], we write
%   norm(m)=sqrt((12-4*trace(R))/(1+trace(R))).
%
%
%
%[1] Förstner, Wrobel, 2004, "Mathematical concepts in Photogrammetry", Ch. 2
%    of "Manual of Photogrammetry", McGlone et al., IAPRS.

t=trace(R);
   
if abs(t-3)<1e-8
        r=zeros(3,1);
        m=zeros(3,1);
        return 
end

nm=sqrt((12-4*t)/(1+t));  %  'nm' is norm of m.
a=-[R(2,3)-R(3,2);R(3,1)-R(1,3);R(1,2)-R(2,1)];
r=a/norm(a) ; % find rotation axes 
m=r*nm;
end


