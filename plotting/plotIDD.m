function plotIDD(varargin)
%PLOTIDD Summary of this function goes here
%   Detailed explanation goes here
figure
hold on

if isstr(varargin{end}) && (strcmp(varargin{end},'RBE') || strcmp(varargin{end},'RBExD'))
    for i = 1:numel(varargin)-1
        plot(matRad_calcIDD(varargin{i}.RBExD),'LineWidth',1.5)
    end
    ylabel('RBExD')
else
    for i = 1:numel(varargin)
        plot(matRad_calcIDD(varargin{i}.physicalDose),'LineWidth',1.5)
    end
    ylabel('physicalDose')
end
xlim([0 250])
end

