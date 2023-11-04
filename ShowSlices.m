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
% Script ShowSlices displays three user-selected brain reconstruction slices.  
% Each slice has two two reconstructions: one from noiseless data, and the 
% other from noisy data.
%
% This script assumes the correct variables in the workspace:
%     X is a set of brain slices listed as vertical vectors in an array.
%     The first 181 are noiseless and the second 181 are noisy.
%
% Author: Tom Asaki
% Version: November 25, 2018
%
%


% choose three slices to display
slices=[50 90 130];

pos=get(0,'screensize');
fw=pos(4)*0.85;
figure('position',[1 1 2*fw/3 fw]);

for k=1:3

  % plot noiseless reconstruction in left figure window
  subplot('position',[0 .333*(k-1) .5 .333])
  imagesc(reshape(X(:,slices(k)),108,108));
  colormap(gray(256));caxis([0 255]);
  set(gca,'Visible','off');

  % plot noisy reconstruction in right figure window
  subplot('position',[.5 .333*(k-1) .5 .333])
  imagesc(reshape(X(:,181+slices(k)),108,108));
  colormap(gray(256));caxis([0 255]);
  set(gca,'Visible','off');

end

