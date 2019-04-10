time_lim_s = 2*60*60;
files = dir('models/*.mat');
for file = files'
    sample(fullfile('models',file.name),time_lim_s);
end