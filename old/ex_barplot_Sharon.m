data_vec=rand(1,20)*0.5;
data_vec2=rand(1,20)*0.5+0.2;

figure;
% You have to plot bar by bar, first the first bar
x_position_bar=1; % position of the bar on the x-axis
data=data_vec;
colorBar='r'; % color of the bar (can be a string or a 3-digit RBG code)
colorError='k';  % color of the error-bar (can be a string or a 3-digit RBG code)
widthBar=0.9; % width of the bar (in units of the x-axis). The bar will be centered in x_position_bar and goes from x_position_bar-widthBar/2 to x_position_bar+widthBar/2
sigFlag=[]; % to perform stats (see comments in the function to use)
widthLine=4; % width of the line of the bar
handle_bar(1)=simpleBarPlot(x_position_bar,data,colorBar,widthBar,colorError,sigFlag,widthLine);

% then second bar
x_position_bar=2; % position of the bar on the x-axis
data=data_vec2;
colorBar=[0.2 1 0.2]; % color of the bar (can be a string or a 3-digit RBG code)
colorError='k';  % color of the error-bar (can be a string or a 3-digit RBG code)
widthBar=0.9; % width of the bar (in units of the x-axis). The bar will be centered in x_position_bar and goes from x_position_bar-widthBar/2 to x_position_bar+widthBar/2
sigFlag=[]; % to perform stats (see comments in the function to use)
widthLine=4; % width of the line of the bar
handle_bar(2)=simpleBarPlot(x_position_bar,data,colorBar,widthBar,colorError,sigFlag,widthLine);

% label the x-axis
xlim([0.2 2.8])
set(gca,'XTick',1:2,'XTickLabel',{'Data 1','Data 2'});
legend(handle_bar,{'Data 1','Data 2'})
ylabel('Whatever dimmension (units)')
format_fig;
title('Example of Bar plot')