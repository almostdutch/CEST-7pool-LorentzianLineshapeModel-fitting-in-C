%
% Program to convert output data from dataFitting from C-style (row major
% order) to MATLAB-style (column major order)
%
% Input:
% outputFileName.mat
%
% Example:
% c2matlab_data_conversion('outputFileName.mat')
%
% Output:
% outputFileName_c2m.mat
%

%
% (c) (2019), Vitaliy Khlebnikov, PhD
%

function c2matlab_data_conversion(InputName)

load(InputName) 

[dim1,dim2,dim3,dim4]=size(fitPars);
Nparameters=size(fitPars,4);
NoffsetsOut=size(fitWaterPool,4);

fitPars=reshape(fitPars,Nparameters,dim3,dim2,dim1);
fitPars=permute(fitPars,numel(size(fitPars)):-1:1);

fitMtPool=reshape(fitMtPool,NoffsetsOut,dim3,dim2,dim1);
fitMtPool=permute(fitMtPool,numel(size(fitMtPool)):-1:1);

fitNoe1Pool=reshape(fitNoe1Pool,NoffsetsOut,dim3,dim2,dim1);
fitNoe1Pool=permute(fitNoe1Pool,numel(size(fitNoe1Pool)):-1:1);

fitNoe2Pool=reshape(fitNoe2Pool,NoffsetsOut,dim3,dim2,dim1);
fitNoe2Pool=permute(fitNoe2Pool,numel(size(fitNoe2Pool)):-1:1);

fitCrPool=reshape(fitCrPool,NoffsetsOut,dim3,dim2,dim1);
fitCrPool=permute(fitCrPool,numel(size(fitCrPool)):-1:1);

fitAmidePool=reshape(fitAmidePool,NoffsetsOut,dim3,dim2,dim1);
fitAmidePool=permute(fitAmidePool,numel(size(fitAmidePool)):-1:1);

fitWaterPool=reshape(fitWaterPool,NoffsetsOut,dim3,dim2,dim1);
fitWaterPool=permute(fitWaterPool,numel(size(fitWaterPool)):-1:1);

save(strcat(InputName(1:end-4),'_c2m.mat'),'fitPars','fitMtPool','fitNoe1Pool','fitNoe2Pool','fitCrPool','fitAmidePool','fitWaterPool','-v7.3')
