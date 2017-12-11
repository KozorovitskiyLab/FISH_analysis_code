clearvars
close all

%What file do you want to read in
[filename, folder]=uigetfile('*.lif');
fullFileName=fullfile(folder,filename);
data=bfopen(fullFileName);


%Where do you want to save output files?
parent_folder=uigetdir;
mkdir(parent_folder,filename);
folder_name=fullfile(parent_folder,filename);

%Run the main program
    FISHCountGeneric(data,folder_name,filename);
