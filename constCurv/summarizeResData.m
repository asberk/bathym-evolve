datapath = 'resdata/';
W = what(datapath);
fileNames = W.mat; % files still have extension 
%[path names ext] = fileparts(W.mat{k});

fprintf('starting summary loop...\n\n');

tic;
for k = 1:length(fileNames) % should be 13 (at the time of writing)
  matFile = [datapath fileNames{k}];
  display(matFile); % should have file extension
  evoLSPlotSummary(matFile);
  fprintf('completed plot summary in %9.6f sec...\n', toc);
end

fprintf('\nfinished loops...\n\n');
fprintf('exiting script...\n');
exit;
