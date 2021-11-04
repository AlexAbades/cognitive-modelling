function NeckerExercise
clear all;
%----------------------------------------------------------------------
% 02458 Cognitive Modelling - Necker Exercise
% Orignial version by Jouko Lampinen, Helsinki University of Technology
% Modified by Tobias Andersen, Technical University of Denmark
%
% Variables:
%  S   - 3D vertex list of the true structure, as 3xN matrix
%  edg - edgea list as Mx2 matrix, each row contains the indices of
%        the vertices linked by the edge
%  I   - projection of S to the image plane, as 2xN matrix
%  M   - projection matrix to project S to I,
%          I=M*[S; ones(1,size(S,2))];
%        M can be constructed by
%          M=viewmtx(AZ,EL);
%	 where Az is the azimuth and EL elevation angle of the observer.
%        Note: for orthographic projection only 3 columns of M are needed,
%        so that I=M(1:2,1:3)*S.
%        For perspective projection (M=viewmtx(AZ,EL,PHI)) M is 2 by 4 and
%        the homogenous coordinates [S;ones(1,size(S,2))] are used.
%        For generality we use full projection matrix here.
%
% The task is to estimate Shat below given the observed image I

%-- vertices of the true underlying scence
S=[ 0 1 1 0 0 1 1 0 ;0 0 1 1 1 1 0 0;0 0 0 0 1 1 1 1]*2-1;

%-- edges as start point, end point index pairs
edg=[ 1 2; 1 4 ; 1 8 ; 2 3 ; 2 7 ; 3 4 ; 3 6 ; 4 5 ; 5 6 ; 5 8 ; 6 7 ; 7 8];

%-- select the viewing angle
AZ=-32; EL=25;
AZ2=-47; EL2=57;
%-- construct the projection matrix

M=viewmtx(AZ,EL); M=M(1:2,:);
M2 = viewmtx(AZ2,EL2); M2=M2(1:2,:);

%-- compute the 2D image
I=M*[S; ones(1,size(S,2))];
I2=M2*[S; ones(1,size(S,2))];

% Make an initial random guess, Sinit, of the true underlying scene
Sinit=rand(size(S));

% Use an optimization routine to find a better guess minimizing
% the error in form of the negative logarithm of the posterior
options = optimset('MaxFunEvals',1000000,'TolFun',1e-3,'TolX',1e-3);
Shat = fminunc(@NeckerError,Sinit,options);

% Draw the true underlying scene, S, and the best fit, Shat.
figure(1); clf;
plotscene(S,edg,'ro-');
hold on; plotscene(Shat,edg,'bo-'); hold off; axis equal; rot(AZ,EL);
title('3D Scene, best fit');axis off;


    function NegLogPost=NeckerError(Sguess)

        %This is the error function that must return the negative logarithm
        %of the posterior probability
        
        %You can use all the variables from the calling function
        %NeckerExercise here to help you
        
        %-- Here is help code for computing the angles between all edges
        anglist=[];
        for k=1:size(Sguess,2),
            con=[edg(find(edg(:,1)==k),2); edg(find(edg(:,2)==k),1)]';
            anglist=[anglist; [ k con([1 2]); k con([1 3]); k con([2 3])]];
        end

        u1=Sguess(:,anglist(:,2))-Sguess(:,anglist(:,1));
        u2=Sguess(:,anglist(:,3))-Sguess(:,anglist(:,1));

        u1=u1./repmat(sqrt(sum(u1.^2)),3,1);
        u2=u2./repmat(sqrt(sum(u2.^2)),3,1);
        Angles=acos(sum(u1.*u2))*180/pi;
        %%-- end of help code
        
        % Get the anlges from the original cube.
        u4=S(:,anglist(:,2))-S(:,anglist(:,1));
        u5=S(:,anglist(:,3))-S(:,anglist(:,1));

        u4=u4./repmat(sqrt(sum(u4.^2)),3,1);
        u5=u5./repmat(sqrt(sum(u5.^2)),3,1);
        Angles_gu = acos(sum(u4.*u5))*180/pi;
        
                
        % We need to compare the original image I with the projection that
        % we've made.
        
        %%%%%%%%%%
        % WE ARE TRYING TO COMPARE THE PERCEPTION THAT HAS A HUMAN BEING
        % WITH ONLY ONE EYE. ONLY ONE POINT OF VIEW
        %%%%%%%%%%
        
        % Without prior:
        %NegLogPost = sum(sum((I-(M*[Sguess; ones(1,size(S,2))])).^2));
        
        % With prior. The more weight sigma has, the less importance that
        % the prior has. Inverse proportion       
        sigma_p = 200000;
        %NegLogPost= sum(sum((I-(M*[Sguess; ones(1,size(S,2))])).^2)) + ...
         %   (1/sigma_p) * sum((Angles-Angles_gu).^2);
        
        %%%%%%%%%
        % IF WE WANT TO SIMULATE THE PERCEPTION OF A HUMAN BEING WITH TWO
        % EYES WE HAVE TO SIMULATE TWO DIFFERENT POINTS OF VIEW. WE CAN
        % MAKE IT WITH TWO DIFFERENT PERCEPTIONS. 
        % We don't need the prior knowladge.
        %%%%%%%%%
        
        
        NegLogPost = sum(sum((I-(M*[Sguess; ones(1,size(S,2))])).^2))+...
            sum(sum((I2-(M2*[Sguess; ones(1,size(S,2))])).^2));
        
        
       
    end

end