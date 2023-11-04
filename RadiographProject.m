%% Radiograph Project - Brain Slicing and Noise Reduction
% Zafran A. Arif - 11594791

% I did this project alongside some of my friends who are very helpful. We
% had a very productive discussion (Dec 3, 2022) that helped us to finish 
% this project.

% All of the problems will be answered through MATLAB comments and will be
% explained in the PDF document report.

% Credits to Zach, Blake, Hideki, and Grace. These are the people in our
% study group.
%% Radiographs
% This particular section was given by Zach but the other sections were
% originally written by myself.

% loading the radiograph data provided by Dr. Tom Asaki
load Lab5Radiographs.mat 
n = 108;
m = 108;
th = linspace(1,179,120); % 1 deg to 179 deg, 120 different views
ScaleFac = 1;
T = tomomap(n,m,th,ScaleFac);
T = full(T);

%% Slices
X = (T'*T)\(T'*B); % array of radiograph constructions
ShowSlices; % code to create slices of specific radiographs

%% First Question
a = length(th); % number or angles taken in radiograph transformation
M = m*a; % verify M = m*a = 12960
N = n*n; % verify N = 108*108

%% Second Question
P = inv(T'*T)*T'; % left inverse transformation
rankP = rank(P); % rank of P

% By Thm 5.1.21, T is injective since we can find the null(T) which is just
% an empty set.

% rank(P) gives us a number 11664 which is equal to N or the dimension of V.
% Since rank(P) is equal to dim(V), we know that P has 11664 linearly
% independent vectors. Thus P must be surjective.

%% Third Question
ShowSlices; % showing the specific slices of the brain
% There's no such thing as the perfect/precise sensor or tool to collect a
% noiseless data. Even with the best radiography equipment that we have today, 
% we still encounter some noises in our data. 

% Thus, it's impossible to get an exact reconstructions.

%% Fourth Question
% Notice that in the radiographs, the noises are present and they were bad for our
% radiograph. These noises are not just in the noisy radiograph, they also
% appear in the noiseless reconstruction as well. This might happen because
% if we read an object using a sensor, we will catch some unwanted data
% outside of the object, thus noisy radiograph.

% We think that there are several possible sources of data noise. Namely,
% Immage Scattering, Dust (or some random particles), Defective Sensors,
% and several errors from the computations.

%% Fifth Question
MatrixA = T'*T; % matrix of T transpose * T (symmetric, diagonalizable by Thm 7.3.12)
[V,D] = eig(MatrixA); % the eigenvalues V and diagonal D
invMatrixA = inv(MatrixA); % Matrix A inverse
eiginvMatrixA = eig(invMatrixA); % the eigenvalues of the inverse of Matrix A

%% Sixth Question
% Let bNoisy = bPerfect + noise, 
% where bNoisy is the Noisy Radiograph, bPerfect is the Perfect Radiograph 
% and noise is the Noise found in Radiograph.

% Let M = Q*D*Q', then inverse of M is invM = Q*invD*Q'

% invM * bNoisy = invM*(bPerfect + noise)
%               = invM*bPerfect + invM*noise
%               = (Q*invD*Q')*(bPerfect)+(Q*inv(D)*Q')*(noise)

% From class, we know that Q*invD*Q' = x, thus
% invM * bNoisy = (Q*invD*Q')*(bPerfect)+(Q*inv(D)*Q')*(noise)
%               = x + (Q*inv(D)*Q')*(noise)
%
% Let xNoise = (Q*inv(D)*Q')*(noise), hence
% invM*bNoisy = x + xNoise

%% Seventh Question
i = 50;
noiseRad = X(:,i+181) - X(:,i); % 50th noisy rad - 50th noiseless rad
normNoise = norm(noiseRad); % magnitude of noise vector k, kth slice
normClean = norm(X(:,i)); % magnitude of clean rad
noiseMatters = normClean - normNoise; % magnitude of noise and clean radiography difference

% Slices: 45, 90, 135
% Noise in the Radiograph: 11664x1 double data
% Noise Magnitude: 6.2522e+03
% Clean Radiograph Magnitude: 8.7171e+03

% The difference (Clean Rad. Magnitude - Noise Magnitude):  2.4648e+03

%% Eighth Question
% Let pseudo-inverse P = sum(VkUk'/sigmak)
% Let Pb = 1/sigma1 * <u1,b>v1 + ... + 1/sigmar<ur,b>*vr

% Suppose the initial radiograph transformation, multiply both sides by P
% P(bNoisy)  = P(bPerfect+noise)
%            = P_bPerfect + P_noise
%            = x + (1/sigma1 * <u1,eta>v1 + ... + 1/sigmar<ur,eta>*vr)

% Notice that sigma_r is a number that is very close to zero. 
% So, if we divide 1 by $\sigma_r$, we will get a very large number. 
% This condition explained the large appearance amount of the noise in our radiograph. 
% In other words, the noises are amplified in the brain reconstruction process. 
% In order to reduce this amplification, we suggest to use a correction. 
% However, the correction method will not give us the perfect answer in the first trial, 
% hence we need to do correction multiple times so that we get the desired clean radiograph. 
% 
% By correction method, we meant enhancing the picture using the nullspace correction method 
% which we found very accurate (we discussed this in class on Friday, December 9). 
% Each iteration of correction will give us a better image result, 
% then after several corrections, we will get an image that is close to the actual object.