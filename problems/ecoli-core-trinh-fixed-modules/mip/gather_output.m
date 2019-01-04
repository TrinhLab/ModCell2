clear;clc;
problem_mip_path = pwd;
modcell_path = fileparts(which('initModCell2.m'));
pin = load(fullfile(modcell_path, 'problems', 'ecoli-core-trinh-fixed-modules','prodnet.mat'));
prodnet = pin.prodnet;

%%
T1 = format_output(fullfile(problem_mip_path,'mc1','mip-lsgcp.csv'),prodnet, 'sGCP');
T2 = format_output(fullfile(problem_mip_path,'mc2','mip-lsgcp.csv'),prodnet, 'sGCP');
T3 = format_output(fullfile(problem_mip_path,'mc3','mip-lsgcp.csv'),prodnet, 'sGCP');

writetable(T1, 'results.xlsx','Sheet','mc1')
writetable(T2, 'results.xlsx','Sheet','mc2')
writetable(T3, 'results.xlsx','Sheet','mc3')
