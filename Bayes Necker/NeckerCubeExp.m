clear all; clc;
%% Infer a 3-D structure of wireframe from a 2-D projection Image	 
%{
S --> Structure. Underlying 3D structure or scene that we can percept from
the 2-D Image (I)
I --> Image in 2D with different interpretations
Bayesian approach to maximaize P(S|I)
For practical reasons we minimize the negative logarithm of the posterior
probability	instead.
M - projection matrix to project S to I, computes a 4-by-4 orthographic or
perspective transformation matrix that projects four-dimensional 
homogeneous vectors onto a two-dimensional view surface (e.g., your 
computer screen).
It's how the prespective of a given object, in this case a cube, changes
due to the observer position. Specifically with the 
%}

%% Necker's Cube 
x = [0  1  1  0  0  0  1  1  0  0  1  1  1  1  0  0];
y = [0  0  1  1  0  0  0  1  1  0  0  0  1  1  1  1];
z = [0  0  0  0  0  1  1  1  1  1  1  0  0  1  1  0];


A = viewmtx(-32,25);
[m,n] = size(x);
x4d = [x(:),y(:),z(:),ones(m*n,1)]';
x2d = A*x4d;
x2 = zeros(m,n); y2 = zeros(m,n);
x2(:) = x2d(1,:);
y2(:) = x2d(2,:);
plot(x2,y2)
axis square


%% 

NeckerExercise