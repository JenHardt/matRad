function plotIDD(varargin)
%PLOTIDD Summary of this function goes here
%   Detailed explanation goes here
figure
hold on

profile = 0;

counter = 0;
for j = 1:length(varargin)
    if isstruct(varargin{j})
        counter = counter + 1;
    end
end

if isstr(varargin{end}) && (strcmp(varargin{end},'RBE') || strcmp(varargin{end},'RBExD'))
    for i = 1:counter
        plot(matRad_calcIDD(varargin{i}.RBExD_MCN,profile),'LineWidth',1.5, 'DisplayName', inputname(i),'LineStyle',varargin{counter+1}{i},'Color',varargin{counter+2}{i})
    end
if profile
    ylabel('profile / RBExD')
 else   
    ylabel('IDD / RBExD')
    end
else
    for i = 1:counter
        plot(matRad_calcIDD(varargin{i}.physicalDose,profile),'LineWidth',1.5, 'DisplayName', inputname(i),'LineStyle',varargin{counter+1}{i},'Color',[varargin{counter+2}{i}])
    end
    if profile
    ylabel('profile / physicalDose')
 else   
    ylabel('IDD / physicalDose')
    end
end
xlabel('depth [voxel]')

% xlim([0 size(varargin{i}.physicalDose,1)*0.8])
legend('Interpreter', 'none','Location','northwest')
end

