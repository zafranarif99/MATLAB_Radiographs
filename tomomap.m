function A=tomomap(n,m,th,ScaleFac)

%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % This file is part of the IMAGE_Math Project and is intended for      %
% % educational use by undergraduate instructors and students. This work %
% % is not for any other use, quotation, or distribution without written %
% % consent of the authors. The authors would like to track usage of     %
% % this code.  Please help us by contacting via info@imagemath.org.     %
% % For the most current version of this code, more information, or      %
% % questions/assistance using the code see http://www.imagemath.org/.   %
% % Copyright 2015-2019.                                                 %    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function TOMOMAP creates a density projection matrix representing the 
% linear transformation from a 2-d square-gridded object space to a series 
% of 1-d projections on a uniform interval spacing.  This is the (linearized) 
% geometric radiography projection.  
%
% Author: Tom Asaki
% Version: July 9, 2018
%
% USAGE:
%
%   A=tomomap(n,m,th)
%   A=tomomap(n,m,th,ScaleFac,SmVal)
%   
% INPUTS:
%
%   n   integer number of cells along one side of the object space.  
%       Objects are defined by values on this n by n grid.
%
%   m   integer number of intervals along the projection plane.
%       This is the number of radiographic pixels per view.
%
%     *** NOTE: The geometric arrays defined by n and m are assumed
%     *** to be co-centered so that the center of the object grid 
%     *** always projects to the center of the radiograph view grid, 
%     *** for each rotation angle.  See drawing below.
%
%   th  ordered vector list of object rotation angles in degrees
%
%   ScaleFac   scale factor for projection intervals.  Object pixels 
%              are fixed at one unit by one unit.  projection intervals 
%              are ScaleFac units.  Default value = 1.
%
%   SmVal      Entries in the computed projection matrix that have
%              absolute value smaller than SmVal are set to zero.
%              Default = 0.001;
% 
% OUTPUTS:
%
%   A   Projection matrix.  A is (m times length(th)) by (n^2).
%       A(i,j) is the fraction of object voxel j that projects to 
%       radiograph pixel i along the perpendicular beam path to pixel i.
%       A is sparse for typical n and m, and is returned in sparse 
%       format.
%
% NOTES:
%
% Example geometric scenario.  In this case, there are two
% views, one at zero degrees and one at 90 degrees.  The object
% array is numbered in Octave/Matlab order.  Radiograph degrees
% are measured east of south (relative to the center of the 
% object array) with individual arrays numbered left to right.  
% In this picture N=n*n.
%
%     +---+------------+---+     -
%     | 1 | n+1        |   |     |
%     +---+            +---+     | 2m
%     | 2 | n+2            |     -
%     +---+    .           |     .
%     | .        .         |     .      radiograph pixels
%     | .          .       |     .        at 90 degrees
%     | .  object voxels   |     -
%     +---+            +---+     |
%     | n |            | N |     | m+1
%     +---+---+--------+---+     -
%     
%
%     |----|----|--   |----|
%       1    2   ...    m
%
%        radiograph pixels 
%         at zero degrees
%
%
%

%%%%% Input Checking and Default Setting %%%%%

A=[];
if    nargin<3 ...
   || n<1 ...
   || m<1 ...
   || round(n)~=n ...
   || round(m)~=m,
   disp('Error: n and m must be positive integers.');
   return
end

ScaleFac_df=1;
SmVal_df=0.001;
if ~exist('ScaleFac','var') || isempty(ScaleFac), ScaleFac=ScaleFac_df; end
if ~exist('SmVal','var') || isempty(SmVal), SmVal=SmVal_df; end


%%%%% Preliminaries %%%%%

% constants
nv=length(th);     % number of views
M=m*nv;            % number of data values
N=n*n;             % number of object locations
pi=4*atan(1);      % pi

% rectify angle set
th=mod(th,360);

% size the operator matrix A
A=sparse(M,N);

% z-coordinates of data 
z=ScaleFac*(-m/2:1:m/2);

% x-coordinates of object as array
xl=(-n/2:1:n/2);
x=repmat(xl,n+1,1);

% y-coordinates of object as array
y=repmat(xl',1,n+1);

%%%%% Main Routine %%%%%

% loop over rotation angles and compute block-rows of A

for k=1:nv

   ang=th(k);

   if ang==0

      BA=zeros(m+1,N);
      for j=1:length(z)
         Q=max(0,x-z(j));
         BA(j,:)=reshape(Q(1:end-1,2:end)-Q(1:end-1,1:end-1),1,[]);
      end
      BA=sparse(BA(1:end-1,:)-BA(2:end,:));

   elseif ang<90

      BA=zeros(m+1,N);
      rad=ang*pi/180;
      factor=(cot(rad)+tan(rad))/2;
      R=[cos(rad) -sin(rad) ; sin(rad) cos(rad)];
      RC=R*[x(:)';y(:)'];
      Rx=reshape(RC(1,:),n+1,n+1);
         
      for j=1:length(z)
         h=Rx-z(j);
         Q=factor*max(0,h).^2;
         BA(j,:)=reshape(Q(1:end-1,2:end)-Q(2:end,2:end)...
                     -Q(1:end-1,1:end-1)+Q(2:end,1:end-1),1,[]);
      end
      BA=sparse(BA(1:end-1,:)-BA(2:end,:));

   elseif ang<180
 
      BA=tomomap(n,m,ang-90,ScaleFac);
      idx=repmat((N-n:-n:0)',1,n)+repmat(1:n,n,1);
      BA=BA(:,idx(:)');

   elseif ang<270

      BA=tomomap(n,m,ang-180,ScaleFac);
      BA=flipud(BA);

   else

      BA=tomomap(n,m,ang-270,ScaleFac);
      idx=repmat((0:n:N-n)',1,n)+repmat(1:n,n,1);
      BA=BA(:,idx(:)');

   end

   A((k-1)*m+1:k*m,:)=BA;
   
end

if SmVal>0
   A(abs(A)<SmVal)=0;
end

return

