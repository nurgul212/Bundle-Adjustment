X = [ -6  -6  -7   0   7   6   6  -3  -3   0   0  -6
      -7   2   1   8   1   2  -7  -7  -2  -2  -7  -7 ;zeros(1,12)];

plot3(X(1,:),X(2,:),X(3,:),'x-')
axis equal

angles=[10,30,230]';
fixed=false;
seq='ats';
v0=[0,0,pi]';

for i=1:angles(1)
    v=[i,0,0]'*pi/180;
    R=EulerRotMat(v0+v,seq,fixed);
    Y=R*X;
    plot3(Y(1,:),Y(2,:),Y(3,:),'x-')
    axis equal
    pause(0.25)
end

for i=1:angles(2)
    v=[angles(1),i,0]'*pi/180;
    R=EulerRotMat(v0+v,seq,fixed);
    Y=R*X;
    plot3(Y(1,:),Y(2,:),Y(3,:),'x-')
    axis equal
    pause(0.25)
end

for i=1:angles(3)
    v=[angles(1),angles(2),i]'*pi/180;
    R=EulerRotMat(v0+v,seq,fixed);
    Y=R*X;
    plot3(Y(1,:),Y(2,:),Y(3,:),'x-')
    axis equal
    pause(0.25)
end
