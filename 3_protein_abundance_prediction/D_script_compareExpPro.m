% Goal:
%   - To compare the exp. vs the predictions
%% Define the save folder
pathway = pwd;
save_dir = 'D';
subfolder = [pathway '\' save_dir];
if ~exist(subfolder, 'dir')
    mkdir(subfolder)
end

%% Load inputs
file = strcat(pathway,'\C\','outerJoin_b_v2'); % v2 is for the whole fly only
plotSave = 1;

%% Run
%% Original
[s_out,data_cell,correlation_coefficient1,rsq,slope,y_int] = graph_scatter_protein_ppm_v3(file,plotSave,subfolder);
%% Move files
s.out = s_out; s.data_cell = data_cell; s.cor = correlation_coefficient1; s.rsq = rsq;
s.slope = slope; s.y_int = y_int;
save(strcat(subfolder,'\','s.mat'),'s');

%% Run
%% log scale
[s_out,data_cell,correlation_coefficient1,rsq,slope,y_int] = graph_scatter_log_protein_ppm(file,plotSave,subfolder);
%% Move files
s_log.out = s_out; s_log.data_cell = data_cell; s_log.cor = correlation_coefficient1; s_log.rsq = rsq;
s_log.slope = slope; s_log.y_int = y_int;
save(strcat(subfolder,'\','s_log.mat'),'s_log');

function [s_out,data_cell,correlation_coefficient1,rsq,...
            slope,y_int] = graph_scatter_log_protein_ppm(file,plotSave,subfolder)

%% Change the param here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ylab = ['log_{10} predicted protein expression (PPM)'];
xlab = ['log_{10} pax-db experimental protein expression (PPM)'];

%% Set the save folder
close all

filename = strcat(file,'.xlsx');
data = readtable(filename);
data_cell = table2cell(data);

x = str2double(string(data_cell(:,2))); %pax-db
y = str2double(string(data_cell(:,end))); % prediction

idx = find(x ~= 0);
x = x(idx); y = y(idx);
idx = find(y ~= 0);
x = x(idx); y = y(idx);
% remove is nan
idx = find(~isnan(x));
x = x(idx); y = y(idx);
% remove is nan
idx = find(~isnan(y));
x = x(idx); y = y(idx);

linearReg = fitlm(x,y);
linearReg.Coefficients
linearReg

s_out = [x,y];
%% Log-log transformation
x2 = log10(x);
y2 = log10(y);

x_bounds = linspace(0,8)';
%% Plot v1 - Original - no Log transformation
linearReg = fitlm(x2,y2);
linearReg.Coefficients
linearReg

%% Pearson test
[R1,P1] = corrcoef(x2,y2,'rows','complete');
correlation_coefficient1 = R1(1,2);
p_value1 = P1(1,2);

[RHO,PVAL] = corr(x2,y2,'Type','Spearman');

%% AVG, SEM
[b_poly,S]= polyfit(x2,y2,1); % 1 is linear. b1 is slope, b2 is y-intercept, S is error estimation structure. This output is primarily used as an input to the polyval function to obtain error estimates.
% Estimate the 95% prediction intervals by using polyconf.
alpha = 0.05; % Significance level
[yfit,delta] = polyval(b_poly,x2,S); % this predict y values. polyval saves typing yourself

slope = b_poly(1); % the first element of b_poly is the slope
y_int = b_poly(2);
x_line= linspace(-4,6);

yresid = y2-yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y2)-1)*var(y2);
rsq = 1-SSresid/SStotal; % This demonstrates that the linear equation b1*x+b2 predicts 53% of the variance in the variable y
     
%% Make 1:1 line
x_lin = log10(linspace(10^-2,10^6,1000));
y_lin = x_lin;

%% Make 1:1 line
% determine how many points are within the bounds
total = length(x2);
diff = sum(abs(log10(y2)-log10(x2)) <= 1); %difference is less than O(1)
percentage = diff./total*100;

x_lin = linspace(-2,6,1000);
y_lin = x_lin;

y_up = x_lin+1;
y_up = x_lin+1;

y_down = x_lin-1;
y_down = x_lin-1;

%% Plot
figure(33)
h1 = scatter(x2,y2,30,'k','X');
hold on
h2 = plot(x_lin,y_lin,'k');
hold on
h3 = plot(x_lin,y_up,'k--');
hold on
h4 = plot(x_lin,y_down,'k--');
hold on
title('Whole adult fly');

grid on
xlabel(xlab)
ylabel(ylab)
xlim([-2 6])
ylim([-2 6])
xticks([-2:2:6])
yticks([-2:2:6])
text(-1.8,5.2,[num2str(diff) ' / ' num2str(total) ' = ' num2str(round(percentage)) '%'],'VerticalAlignment','top', 'HorizontalAlignment','left',...
    'FontSize',12);
set(gca,'FontSize',13)

if plotSave == 1
    filename_image_tiff = strcat(subfolder, '\v2_log.tiff');
    print([filename_image_tiff],'-dtiffn','-r300')% saveas(gcf,'shortVSlong.tiff');
    saveas(gcf,strcat(subfolder, '\v2_log.fig'))% saveas(gcf,'shortVSlong.tiff');
    exportgraphics(gcf,strcat(subfolder, '\v2_log.pdf'))% saveas(gcf,'shortVSlong.tiff');

end

end 


function [s_out,data_cell,correlation_coefficient1,rsq,...
            slope,y_int] = graph_scatter_protein_ppm_v3(file,plotSave,subfolder)
%% Change the param here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlab = ['Protein abundance, experimental (PPM)'];
ylab = ['Protein abundance, predicted (PPM)'];

%% Set the save folder
close all

filename = strcat(file,'.xlsx');
data = readtable(filename);
data_cell = table2cell(data);

x = str2double(string(data_cell(:,2))); %pax-db
y = str2double(string(data_cell(:,end))); %  prediction

idx = find(x ~= 0);
x = x(idx); y = y(idx);
idx = find(y ~= 0);
x = x(idx); y = y(idx);
% remove is nan
idx = find(~isnan(y));
x = x(idx); y = y(idx);

linearReg = fitlm(x,y);
linearReg.Coefficients
linearReg

s_out = [x,y];
%% Log-log transformation
x2 = x;
y2 = y;

x_bounds = linspace(10^-2,10^8)';
%% Plot v1 - Original - no Log transformation
linearReg = fitlm(x2,y2);
linearReg.Coefficients
linearReg

%% Pearson test
[R1,P1] = corrcoef(x2,y2,'rows','complete');
correlation_coefficient1 = R1(1,2);
p_value1 = P1(1,2);
[RHO,PVAL] = corr(x2,y2,'Type','Spearman');

%% AVG, SEM
[b_poly,S]= polyfit(x2,y2,1); % 1 is linear. b1 is slope, b2 is y-intercept, S is error estimation structure. This output is primarily used as an input to the polyval function to obtain error estimates.
% Estimate the 95% prediction intervals by using polyconf.
alpha = 0.05; % Significance level
[yfit,delta] = polyconf(b_poly,x2,S,'alpha',alpha);
[yfit,delta] = polyval(b_poly,x2,S); % this predict y values. polyval saves typing yourself
slope = b_poly(1); % the first element of b_poly is the slope
y_int = b_poly(2);
x_line= linspace(-4,6);
yresid = y2-yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y2)-1)*var(y2);
rsq = 1-SSresid/SStotal; % This demonstrates that the linear equation b1*x+b2 predicts 53% of the variance in the variable y
%% Make 1:1 line
% determine how many points are within the bounds
total = length(x2);
diff = sum(abs(log10(y2)-log10(x2)) <= 1); %difference is less than O(1)
percentage = diff./total*100;
x_lin = linspace(10^-2,10^6,1000);
y_lin = x_lin;
y_up = 10.*x_lin;
y_down = x_lin./10;

%% 2023. 11. 15. 
% Check how many values are within the 95% confidence interval
y2_upper_2delta = yfit+2*delta;
y2_lower_2delta = yfit-2*delta;
within_interval = (y2 >= y2_lower_2delta) & (y2 <=y2_upper_2delta);
within_interval_sum = sum(within_interval); 
% Count the number of values within the interval
percentage_within_interval = (within_interval_sum / length(x2)) * 100;
disp([newline 'Number of values within the 95% confidence interval: ' num2str(within_interval_sum),...
            ' out of ' num2str(length(x2)) ' (' num2str(percentage_within_interval) ' %)']);
%% Plot
figure(33)
h1 = scatter(x2,y2,30,'k','X');
grid on
hold on
h2 = plot(x_lin,y_lin,'k');
hold on
h3 = plot(x_lin,y_up,'k--');
hold on
h4 = plot(x_lin,y_down,'k--');
hold on
title('Whole adult fly');
set(gca,'xscale','log','yscale','log');
xlabel(xlab)
ylabel(ylab)
xlim([10^-2 10^6])
ylim([10^-2 10^6])
text(10^-1.8,10^5.8,[num2str(diff) ' / ' num2str(total) ' = ' num2str(round(percentage)) '%'],'VerticalAlignment','top', 'HorizontalAlignment','left',...
    'FontSize',12);
set(gca,'FontSize',13)
grid on
if plotSave == 1
    filename_image_tiff = strcat(subfolder, '\v2.tiff');
    print([filename_image_tiff],'-dtiffn','-r300')% saveas(gcf,'shortVSlong.tiff');
    exportgraphics(gcf,strcat(subfolder, '\v2.pdf'))% saveas(gcf,'shortVSlong.tiff');
    saveas(gcf,strcat(subfolder, '\v2.fig'))% saveas(gcf,'shortVSlong.tiff');

end
end 
