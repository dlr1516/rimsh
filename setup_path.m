function setup_path()
%SETUP_PATH  Add the library folders of this repository to the MATLAB path.
%
%   Run once per session from the repository root:
%       setup_path
%   Use SAVEPATH afterwards to make the change permanent.

  root = fileparts(mfilename('fullpath'));
  addpath(genpath(fullfile(root, 'src')));
  fprintf('Added %s (and subfolders) to the MATLAB path.\n', fullfile(root, 'src'));
end
