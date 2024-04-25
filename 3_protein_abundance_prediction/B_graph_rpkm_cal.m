%% Define the save folder
pathway = pwd;
save_dir = 'B';
subfolder = [pathway '\' save_dir];
if ~exist(subfolder, 'dir')
    mkdir(subfolder)
end

%% Load inputs
file = strcat(pathway,'\A\','results_gene2param_T2uuu'); % v2 is for the whole fly only
plotSave = 1;

%% Run
[s_out_log,data_cell_log,correlation_coefficient1_log,slope,y_int] = graph_scatter_log(file,plotSave,subfolder);
s_log.out = s_out_log; s_log.data_cell_log = data_cell_log;
s_log.data_cell_log = correlation_coefficient1_log;
save(strcat(subfolder,'\','s_log.mat'),'s_log');

[s_out,data_cell,correlation_coefficient1,slope,y_int] = graph_scatter(file,plotSave,subfolder)
s.out = s_out; s.data_cell = data_cell; s.cor = correlation_coefficient1;
s.data_cell = data_cell; s.y_int = y_int;
s.slope = slope;
save(strcat(subfolder,'\','s.mat'),'s');

function [s_out,data_cell,correlation_coefficient1,slope,y_int] = graph_scatter_log(file,plotSave,subfolder)
% This function will graph the basic scatter plots with the R^2 values
% input:
%   - 2 columns => [x, y]
xlab = 'Gene expression (log_{10}RPKM)';
ylab = 'Protein abundance, predicted (log_{10}PPM)';

%% Set the save folder
close all
filename = strcat(file,'.xlsx');
data = readtable(filename);
data_cell = table2cell(data);

x = str2double(string(data_cell(:,end-2)));
y = str2double(string(data_cell(:,end)));

idx = find(x ~= 0);
x = x(idx); y = y(idx);
idx = find(y ~= 0);
x = x(idx); y = y(idx);
% remove is nan
idx = find(~isnan(x));
x = x(idx); y = y(idx);

total_nucleotides = 2.6e-15;
avo = 6.022e23;
rpkm = x./total_nucleotides./avo.*10^6.*10^3;
x = rpkm;

linearReg = fitlm(x,y);
linearReg.Coefficients
linearReg

s_out = [x,y];
%% Log-log transformation
x2 = log10(x);
y2 = log10(y);

x_bounds = linspace(0,8)';
% x_bounds = linspace(10^-2,10^8)';

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

[x2_sorted,idx] = sort(x2,'ascend');
yfit_sorted = yfit(idx);
delta_sorted = delta(idx);

%% Confidence interval
slope = b_poly(1); % the first element of b_poly is the slope
y_int = b_poly(2);
x_line= linspace(-4,6);

yresid = y2-yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y2)-1)*var(y2);
rsq = 1-SSresid/SStotal; % This demonstrates that the linear equation b1*x+b2 predicts 53% of the variance in the variable y

%% Make 1:1 line
x_lin = log10(linspace(10^-2,10^6,1000))';
% y_lin = x_lin;
y_lin = slope.*x_lin+y_int;
[y_lin_fit,delta_lin_fit] = polyconf(b_poly,x_lin,S,'alpha',alpha);

%% Plot
figure(33)
h1 = scatter(x2,y2,30,'k','X');
hold on
xlabel(xlab)
ylabel(ylab)
xlim([-1 4])
ylim([-2 6])
xticks([-2:2:6])
yticks([-2:2:6])
text(-0.8,5.8,['\rho = ' num2str(round(correlation_coefficient1,2))], 'VerticalAlignment','top', 'HorizontalAlignment','left',...
    'FontSize',11)
set(gca,'FontSize',13)
grid on

if plotSave == 1
    filename_image_tiff = strcat(subfolder, '\v2.tiff');
    print([filename_image_tiff],'-dtiffn','-r300')% saveas(gcf,'shortVSlong.tiff');
end
end 



function [s_out,data_cell,correlation_coefficient1,slope,y_int] = graph_scatter(file,plotSave,subfolder)
xlab = 'Gene expression (RPKM)';
ylab = ['Protein abundance, predicted (PPM)'];

%% Set the save folder
close all

filename = strcat(file,'.xlsx');
data = readtable(filename);
data_cell = table2cell(data);

x = str2double(string(data_cell(:,end-2)));
y = str2double(string(data_cell(:,end)));

idx = find(x ~= 0);
x = x(idx); y = y(idx);
idx = find(y ~= 0);
x = x(idx); y = y(idx);
% remove is nan
idx = find(~isnan(x));
x = x(idx); y = y(idx);

total_nucleotides = 2.6e-15;
avo = 6.022e23;
rpkm = x./total_nucleotides./avo.*10^6.*10^3;
x = rpkm;

linearReg = fitlm(x,y);
linearReg.Coefficients
linearReg

s_out = [x,y];

x2 = (x);
y2 = (y);

x_bounds = linspace(0,8)';

%% Plot v1 - Original - no Log transformation
linearReg = fitlm(x2,y2);
linearReg.Coefficients
linearReg

%% Make 1:1 line
% determine how many points are within the bounds
total = length(x2);
diff = sum(abs(log10(y2)-log10(x2)) <= 1); %difference is less than O(1)
percentage = diff./total*100;
x_lin = linspace(10^-2,10^6,1000);
y_lin = x_lin;
y_up = 10.*x_lin;
y_down = x_lin./10;

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
     
%% Plot
figure(33)
h1 = scatter(x2,y2,30,'k','X');
hold on
grid on
xlabel(xlab)
ylabel(ylab)
xlim([10^-1 10^4])
ylim([10^-2 10^6])
text(10^-0.8,10^5.8,['\rho = ' num2str(round(correlation_coefficient1,2))], 'VerticalAlignment','top', 'HorizontalAlignment','left',...
    'FontSize',12)
set(gca,'FontSize',13)
set(gca,'xscale','log')
set(gca,'yscale','log')
% hold off
if plotSave == 1
    filename_image_tiff = strcat(subfolder, '\v2_noLogRPKM.tiff');
    print([filename_image_tiff],'-dtiffn','-r300')% saveas(gcf,'shortVSlong.tiff');
end

end 
