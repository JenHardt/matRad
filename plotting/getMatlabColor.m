function rgb = getMatlabColor(color)

if strcmp(color,'blue')
    rgb = [0, 0.4470, 0.7410];
elseif strcmp(color,'orange')
    rgb = [0.8500, 0.3250, 0.0980];
elseif strcmp(color,'yellow')
    rgb = [0.9290, 0.6940, 0.1250];
elseif strcmp(color,'purple')
    rgb = [0.4940, 0.1840, 0.5560];
elseif strcmp(color,'green')
    rgb = [0.4660, 0.6740, 0.1880];
elseif strcmp(color,'cyan')
    rgb = [0.3010, 0.7450, 0.9330];
elseif strcmp(color,'red')
    rgb = [0.6350, 0.0780, 0.1840];
end

end

