% Start with a folder and get a list of all subfolders.
% Finds and prints names of all files in 
% that folder and all of its subfolders.
% Similar to imageSet() function in the Computer Vision System Toolbox: http://www.mathworks.com/help/vision/ref/imageset-class.html
clc;    % Clear the command window.
workspace;  % Make sure the workspace panel is showing.
format long g;
format compact;

% Define a starting folder.
start_path = 'E:\MATLAB\TDTbin2mat\Inscopix Python\Individual_Mouse_Cell_Data'; %fullfile(matlabroot, '\toolbox');
if ~exist(start_path, 'dir')
	start_path = matlabroot;
end
% Ask user to confirm or change.
uiwait(msgbox('Pick a starting folder on the next window that will come up.'));
topLevelFolder = uigetdir(start_path);
if topLevelFolder == 0
	return;
end
% Get list of all subfolders.
allSubFolders = genpath(topLevelFolder);
% Parse into a cell array.
remain = allSubFolders;
% Initialize variable that will contain the list of filenames so that we can concatenate to it.
listOfFolderNames = {};
while true
	[singleSubFolder, remain] = strtok(remain, ';');
	if isempty(singleSubFolder)
		break;
	end
	listOfFolderNames = [listOfFolderNames singleSubFolder];
end
numberOfFolders = length(listOfFolderNames)

% Process all image files in those folders.
for k = 1 : numberOfFolders
	% Get this folder and print it out.
	thisFolder = listOfFolderNames{k};
	fprintf('Processing folder %s\n', thisFolder);
	
	% Get ALL files.
	filePattern = sprintf('%s/*.*', thisFolder);
	baseFileNames = dir(filePattern);
% 	% Get m files.
% 	filePattern = sprintf('%s/*.m', thisFolder);
% 	baseFileNames = dir(filePattern);
% 	% Add on FIG files.
% 	filePattern = sprintf('%s/*.fig', thisFolder);
% 	baseFileNames = [baseFileNames; dir(filePattern)];
	% Now we have a list of all files in this folder.
	
	numberOfImageFiles = length(baseFileNames);
	if numberOfImageFiles >= 1
		% Go through all those files.
		for f = 1 : numberOfImageFiles
			fullFileName = fullfile(thisFolder, baseFileNames(f).name);
			fprintf('     Processing file %s\n', fullFileName);
		end
	else
		fprintf('     Folder %s has no files in it.\n', thisFolder);
	end
end

