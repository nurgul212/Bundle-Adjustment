function [P2,X]=nister5pt(pt1n,pt2n,relax)
%NISTER5PT Estimate relative orientation using Nister's 5-pt algorithm.
%
%   [P2,X]=NISTER5PT(P,Q), where P and Q are 3-by-N arrays with
%   normalized homogeneous points in the left and right camera,
%   respectively, returns the right camera matrix P2 and the estimated 3D
%   points X.
%
%   ...=NISTER5PT(P,Q,TRUE) relaxes the 'no imaginary component'
%   constraint in ESSMAT5.
%
%See also: ESSMAT5, CAMSFROME.

% $Id: nister5pt.m 1287 2014-12-08 15:22:35Z niclas $

if nargin<3, relax=false; end

% Solve for E matrices.
EE=essmat5(pt2n,pt1n,relax);

% Check each possible E matrix.
PP2=zeros(3,4,size(EE,2));
XX=cell(1,size(EE,2));
SP=nan(4,size(EE,2));

for i=1:size(EE,2)
    % Determine best P1/P2 combination from this essential matrix.
    [P1,P2,X,inFront,sp]=camsfrome(reshape(EE(:,i),3,3),pt1n,pt2n,-1);
    PP2(:,:,i)=P2;
    XX{i}=X;
    SP(:,i)=sp;
    if 0
        C2=euclidean(null(P2));
        R2=P2(:,1:3);
        ang=derotmat3d(R2)';
        EEO=[zeros(7,1),[C2;ang;0]];
        fig=figure(tagfigure('nister'));
        T0=eye(4);
        camSize=0.1;
        pm_plotmulti([],EEO,X(1:3,:),false(1,size(X,2)),true(1,2),...
                     [],[],[],[],[],[],[],{},camSize,fig,T0);
        pause
    end
end

% Find best P1/P2 among all.

% Use number of points in front of the cameras as first criterion.
D=SP(1,:)-SP(2,:);
best=find(D==max(D));

if length(best)==1
    P2=PP2(:,:,best);
    X=XX{best};
else
    % Multiple solutions with same # of points in front of cameras.
    % Determine which one has the smallest reprojection error.
    
    bestErr=inf;
    bestP2=nan(3,4);
    bestX=[];

    % Pre-Euclideanize.
    pt1nE=euclidean(pt1n);
    pt2nE=euclidean(pt2n);
    for i=1:length(best)
        ci=best(i);
        P1=eye(3,4);
        P2=PP2(:,:,ci);
        X=homogeneous(XX{ci});
        % Reprojection error.
        err=sum(sum((euclidean(P1*X)-pt1nE).^2))+...
            sum(sum((euclidean(P2*X)-pt2nE).^2));
        if err<bestErr
            bestErr=err;
            bestP2=P2;
            bestX=XX{ci};
        end
    end
    X=bestX;
    P2=bestP2;
end
