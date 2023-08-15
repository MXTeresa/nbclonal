function void = main_FigureS2_plot(void)

clear all; close all; clc;
infile = 'figure_SARSCOV2_T005'; load(infile);
thresh = 1; R0val = 1;
infileS1_IAV = strcat('figureS2_SC2_results_thres', int2str(thresh), '_R0', int2str(R0val)); load(infileS1_IAV);
paneltitle = 'R_0 = 2.6; threshold = 0.5%'
plot_figureS1_panel(params, results, results_emp, 3, 8, 1, 2, paneltitle);

clear all; clc;
infile = 'figure_SARSCOV2_T03'; load(infile);
thresh = 2; R0val = 1;
infileS1_IAV = strcat('figureS2_SC2_results_thres', int2str(thresh), '_R0', int2str(R0val)); load(infileS1_IAV);
paneltitle = 'R_0 = 2.6; threshold = 3%'
plot_figureS1_panel(params, results, results_emp, 3, 8, 4, 5, paneltitle);

clear all; clc;
infile = 'figure_SARSCOV2_T07'; load(infile);
thresh = 3; R0val = 1;
infileS1_IAV = strcat('figureS2_SC2_results_thres', int2str(thresh), '_R0', int2str(R0val)); load(infileS1_IAV);
paneltitle = 'R_0 = 2.6; threshold = 7%'
plot_figureS1_panel(params, results, results_emp, 3, 8, 7, 8, paneltitle);

clear all; clc;
infile = 'figure_SARSCOV2_T005'; load(infile);
thresh = 1; R0val = 2;
infileS1_IAV = strcat('figureS2_SC2_results_thres', int2str(thresh), '_R0', int2str(R0val)); load(infileS1_IAV);
paneltitle = 'R_0 = 7.4; threshold = 0.5%'
plot_figureS1_panel(params, results, results_emp, 3, 8, 9, 10, paneltitle);

clear all; clc;
infile = 'figure_SARSCOV2_T03'; load(infile);
thresh = 2; R0val = 2;
infileS1_IAV = strcat('figureS2_SC2_results_thres', int2str(thresh), '_R0', int2str(R0val)); load(infileS1_IAV);
paneltitle = 'R_0 = 7.4; threshold = 3%'
plot_figureS1_panel(params, results, results_emp, 3, 8, 12, 13, paneltitle);

clear all; clc;
infile = 'figure_SARSCOV2_T07'; load(infile);
thresh = 3; R0val = 2;
infileS1_IAV = strcat('figureS2_SC2_results_thres', int2str(thresh), '_R0', int2str(R0val)); load(infileS1_IAV);
paneltitle = 'R_0 = 7.4; threshold = 7%'
plot_figureS1_panel(params, results, results_emp, 3, 8, 15, 16, paneltitle);

clear all; clc;
infile = 'figure_SARSCOV2_T005'; load(infile);
thresh = 1; R0val = 3;
infileS1_IAV = strcat('figureS2_SC2_results_thres', int2str(thresh), '_R0', int2str(R0val)); load(infileS1_IAV);
paneltitle = 'R_0 = 14.9; threshold = 0.5%'
plot_figureS1_panel(params, results, results_emp, 3, 8, 17, 18, paneltitle);

clear all; clc;
infile = 'figure_SARSCOV2_T03'; load(infile);
thresh = 2; R0val = 3;
infileS1_IAV = strcat('figureS2_SC2_results_thres', int2str(thresh), '_R0', int2str(R0val)); load(infileS1_IAV);
paneltitle = 'R_0 = 14.9; threshold = 3%'
plot_figureS1_panel(params, results, results_emp, 3, 8, 20, 21, paneltitle);

clear all; clc;
infile = 'figure_SARSCOV2_T07'; load(infile);
thresh = 3; R0val = 3;
infileS1_IAV = strcat('figureS2_SC2_results_thres', int2str(thresh), '_R0', int2str(R0val)); load(infileS1_IAV);
paneltitle = 'R_0 = 14.9; threshold = 7%'
plot_figureS1_panel(params, results, results_emp, 3, 8, 23, 24, paneltitle);
