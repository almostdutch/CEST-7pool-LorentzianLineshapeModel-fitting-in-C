function checkData(InputName)
%
% Program to visualize fitted maps
%
% Input:
% outputFileName_c2m.mat
%
% Example:
% checkData('outputFileName_c2m.mat')
%

load(InputName, 'fitPars')

n=[2 5 8 11 14 17];
names={'MT pools' 'NOE1 pool', 'NOE2 pool', 'Cr pool', 'Amide pool', 'Water pool'};

figure
for ii=1:numel(n)
    subplot(2,3,ii)
    imagesc(fitPars(:,:,1,n(ii))),axis off, axis equal, colormap parula
    title(names(ii))
end
end