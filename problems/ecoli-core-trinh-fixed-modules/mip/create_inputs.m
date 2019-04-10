%% All products
modcell_path = fileparts(which('initModCell2.m'));
pin = load(fullfile(modcell_path, 'problems', 'ecoli-core-trinh-fixed-modules','prodnet.mat'));

prodnet2mip('mip-lsgcp',pin.prodnet, 'lsGCP')
