function FISHCountGeneric(data,destination_folder,filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%What stacks do you want to analyze?
stackSelector={'Which stacks do you want to analyze?'};
selectedStacks=inputdlg(stackSelector);
%Str2num is required to parse commas correctly
selectedStacks=str2num(selectedStacks{:});

% Total channels you have imaged?
totalOfChannels=questdlg('How many channels you have imaged?','Channel Number','2','3','4','4');

%which channels you want to analyze?
channelSelector={'Which channels (more than 1) do you want to analyze?'};
selectedchannel=inputdlg(channelSelector);
selectedchannel = str2num(selectedchannel{:}); %#ok<ST2NM>

numberOfChannels = length(selectedchannel);
originalChannelType=cell(numberOfChannels,1);
originalChannelType{numberOfChannels,1}=[];
threshold=zeros(numberOfChannels,1);

channelType=cell(numberOfChannels,1);
channelType{numberOfChannels,1}=[];


%What type is each channel?
for i=1:numberOfChannels
    j=num2str(selectedchannel(i));
    originalChannelType{i,1}=questdlg(['What type of channel is ', j, '?'], 'Channel Type', 'Fill', 'Not Fill', 'Not Fill');
end

fillType=questdlg('Is the fill nuclear or somatic?', 'Fill type', 'Nucleus','Soma','Nucleus');

%Ask user to choose z planes to include
zPlaneNumberChoice={'Enter total number of z planes'};
dlg_title1='Choose how many z planes to include.';
num_lines1=1;
defaultans1={'6'};
zPlaneTotal=str2double(inputdlg(zPlaneNumberChoice,dlg_title1,num_lines1,defaultans1));
%Initialize some variables
splitChannelData=cell(zPlaneTotal,numberOfChannels);
splitChannelData{zPlaneTotal,numberOfChannels}=[];
mask=cell(numberOfChannels,1);
mask{numberOfChannels,1}=[];
maskHold=cell(numberOfChannels,1);
maskHold{numberOfChannels,1}=[];

%stacks
numberOfStacks=length(selectedStacks);

numberOfCategories=(factorial(numberOfChannels-1))+1;
categoryCount=zeros(numberOfStacks,numberOfCategories);



for i=1:numberOfStacks
    k = selectedStacks(i);
    %We have now gotten this down to a single stack for analysis
    totalNumberOfZPlanes=length(data);
    
    %What z-planes do you want to include?
    %Choose a starting plane
    zPlaneStartChoice={'Enter first z plane of stack' num2str(k)};
    zPlaneStartChoice=strjoin(zPlaneStartChoice);
    dlg_title2='Enter the first z plane.';
    num_lines2=1;
    defaultans2={'2'};
    zPlaneStart=str2double(inputdlg(zPlaneStartChoice,dlg_title2,num_lines2,defaultans2));
    
    
    %Assign the images to individual channels
    for l = 1:zPlaneTotal   % firstZPlane:lastZPlane
        for cc = 1:numberOfChannels
            Plane2Channel = (zPlaneStart - 1)* str2double(totalOfChannels) + (l -1)* str2double(totalOfChannels) + selectedchannel(cc);
            splitChannelData{l,cc}=data{k,1}{Plane2Channel,1};
        end
    end
    
    %Create the sums/averages of the different channel types for later
    %mask generation
    FillChannelCount=0;
    
    %I think the code below could probably be condensed dramatically
    for n=1:numberOfChannels
        
        if strcmp(originalChannelType{n,1},'Fill')
            
            channel3D=cat(3,splitChannelData{:,n});
            mask{end}=uint8(mean(channel3D,3));
            %This could be a good time to choose a threshold
            maskHold{end}=mask{end};
            image1=imshow(mask{end});
            image2=imcontrast(gca);
            uiwait(msgbox('Please change the Window minimum value to the appropriate threshold, and click in another field to register the change. Then, click OK.'));
            FillThresholds=get(gca,'CLim');
            %Line below is hard coded for channel 4 DAPI
            threshold(end)=FillThresholds(1);
            
            close all
            mask{end}(mask{end}<threshold(end))=0;
            channelType{end}=originalChannelType{n,1};
            
            FillChannelCount=FillChannelCount+1;
            
            
        elseif strcmp(originalChannelType{n,1},'Not Fill')
            channel3D=cat(3,splitChannelData{:,n});
            mask{n-FillChannelCount}=uint8(mean(channel3D,3));
            
            maskHold{n-FillChannelCount}=mask{n-FillChannelCount};
            image1=imshow(mask{n-FillChannelCount});
            image2=imcontrast(gca);
            uiwait(msgbox('Please change the Window minimum value to the appropriate threshold, and click in another field to register the change. Then, click OK.'));
            lushThresholds=get(gca,'CLim');
            threshold(n-FillChannelCount)=lushThresholds(1);
            close all
            channelType{n-FillChannelCount}=originalChannelType{n,1};
            
            mask{n-FillChannelCount}(mask{n-FillChannelCount}<threshold(n-FillChannelCount))=0;
            
            
        end
        
    end
    
    %Now generate the initial files for the masks
    fullMask=uint8(sum(cat(3,mask{:}),3));
    FillMask=uint8(sum(cat(3,mask{end}),3));
    %There is no need for the FISH mask
    %FISHMask=uint8(sum(cat(3,mask{1:3}),3));
    
    
    if strcmp(fillType,'Nucleus')
        
        %CellSegmenter will do some image processing on the data, and then
        %perform a watershed
        Ld2=CellSegmenterNucleus(fullMask,FillMask,FillThresholds,threshold);
        
    elseif strcmp(fillType,'Soma')
        Ld2=CellSegmenterSoma(fullMask,FillMask,FillThresholds,threshold);
    end
    
    numberOfCells=max(max(Ld2));
    punctaPerCell=zeros(numberOfCells,1);
    allChannelsPunctaPerCell=zeros(numberOfCells,numberOfChannels-1);
    punctaCutoff=zeros(1,numberOfChannels-1);
    logicalPunctaPerCell=zeros(numberOfCells,numberOfChannels-1);
    outputPunctaPerCell = zeros(numberOfCells,numberOfChannels-1);
    for n=1:(numberOfChannels-1)
        %FastPeakFind2 with initialization of threshold set
        peakThreshold=threshold(n);
        [punctaPositions,peakMap]=FastPeakFind5(mask{n},peakThreshold,fspecial('average',2),1,channelType{n},1);
        
        %Punctacount
        logicalPeakMap=logical(peakMap);
        punctaCount=Ld2(logicalPeakMap);
        
        for q=1:numberOfCells
            punctaPerCell(q,1)=sum(punctaCount==q);
            
            testPunctaPerCellHoldAllNew{i,n}(q,1)=punctaPerCell(q);
        end
        
        mostPuncta=max(punctaPerCell);
        edges= (0.5:1:mostPuncta+0.5);
        %Display histogram and ask user for threshold that counts as
        %expression
        figure, hist(punctaPerCell,edges)
        punctaCutoffChoice={'Enter minimum number of puncta that counts as expression.'};
        dlg_title2='Puncta Cutoff';
        num_lines2=1;
        defaultans2={'5'};
        punctaCutoff(1,n)=str2double(inputdlg(punctaCutoffChoice,dlg_title2,num_lines2,defaultans2));
        punctaPerCell(punctaPerCell<punctaCutoff(1,n))=0;
        logicalPunctaPerCell(:,n)=logical(punctaPerCell);
        outputPunctaPerCell(:,n) = testPunctaPerCellHoldAllNew{i,n};
    end
    
    %Now produce a file with a user entered name and directory location
    full_filename1=fullfile(destination_folder, [filename 'Stack ' num2str(k) ' Puncta Per Cell.xlsx']);
    fullworkspace_filename=fullfile(destination_folder,[filename 'Stack ' num2str(k) ' workspace.mat']);
    
    
    %Close all figures before proceeding to save memory
    delete(findall(0,'Type','figure'))
    
    %Code below is not elegant, but is functional
    categoryCount(i,1)=selectedStacks(i);
    for c=1:numberOfCategories
        if numberOfChannels==2
            if c==1
                categoryCount(i,c+1)=sum((logicalPunctaPerCell(:,1)==true));
            elseif c==2
                categoryCount(i,c+1)=sum((logicalPunctaPerCell(:,1)==false));
            end
        else
            if numberOfChannels==3
                if c==1
                    categoryCount(i,c+1)= sum((logicalPunctaPerCell(:,1)==true)&(logicalPunctaPerCell(:,2)==false));
                elseif c==2
                    categoryCount(i,c+1)= sum((logicalPunctaPerCell(:,1)==false)&(logicalPunctaPerCell(:,2)==true));
                elseif c==3
                    categoryCount(i,c+1)= sum((logicalPunctaPerCell(:,1)==true)&(logicalPunctaPerCell(:,2)==true));
                end
            else
                if c==1
                    
                    categoryCount(i,c+1)= sum((logicalPunctaPerCell(:,1)==true)&(logicalPunctaPerCell(:,2)==false)&(logicalPunctaPerCell(:,3)==false));
                elseif c==2
                    categoryCount(i,c+1)= sum((logicalPunctaPerCell(:,1)==false)&(logicalPunctaPerCell(:,2)==true)&(logicalPunctaPerCell(:,3)==false));
                elseif c==3
                    categoryCount(i,c+1)= sum((logicalPunctaPerCell(:,1)==false)&(logicalPunctaPerCell(:,2)==false)&(logicalPunctaPerCell(:,3)==true));
                elseif c==4
                    categoryCount(i,c+1)= sum((logicalPunctaPerCell(:,1)==true)&(logicalPunctaPerCell(:,2)==true)&(logicalPunctaPerCell(:,3)==false));
                elseif c==5
                    categoryCount(i,c+1)= sum((logicalPunctaPerCell(:,1)==true)&(logicalPunctaPerCell(:,2)==false)&(logicalPunctaPerCell(:,3)==true));
                elseif c==6
                    categoryCount(i,c+1)= sum((logicalPunctaPerCell(:,1)==false)&(logicalPunctaPerCell(:,2)==true)&(logicalPunctaPerCell(:,3)==true));
                elseif c==7
                    categoryCount(i,c+1)= sum((logicalPunctaPerCell(:,1)==true)&(logicalPunctaPerCell(:,2)==true)&(logicalPunctaPerCell(:,3)==true));
                end
                
            end
        end
    end
    
    save(fullworkspace_filename, 'channelType','FillMask','FillThresholds','filename','Ld2','logicalPeakMap','logicalPunctaPerCell','mask','numberOfCells','originalChannelType','outputPunctaPerCell','peakMap','punctaCount','punctaCutoff','splitChannelData','stackSelector','testPunctaPerCellHoldAllNew','threshold','zPlaneStart','zPlaneTotal','categoryCount');
    
    
end
%Assign the total data to the workspace. Note: this is not ideal.
assignin('base','cellsByCategory',categoryCount);

%The code below can be uncommented if you want to save as an excel file
%textStacks=num2str(selectedStacks);
%full_filename=fullfile(destination_folder, [filename 'Stacks ' textStacks ' Cell Proportions.xlsx']);
%Now save all the data into it. B2 sets the row and column in excel that
%the data will start at
%xlswrite(full_filename, selectedStacks,1,'A2', categoryCount,1, 'B2');
%xlswrite(full_filename, cellsByCategory,1, 'B2');
%xlswrite(full_filename, categoryCount,1, 'B2');


%The code below can be uncommented if you want to specify where you save your .mat file. 
%clear('data');
%[workspace_filename, path]=uiputfile('*.mat','Save Data As');
%fullworkspace_filename=fullfile(path,workspace_filename);
%save(fullworkspace_filename)
%save(fullworkspace_filename, 'C_data', 'cellsByCategory','channelType','FillMask','DAPIThresholds','filename','Ld2','logicalPeakMap','logicalPunctaPerCell','mask','numberOfCells','originalChannelType','outputPunctaPerCell','peakMap','punctaCount','punctaCutoff','splitChannelData','stackSelector','testPunctaPerCellHoldAllNew','threshold','zPlaneStart','zPlaneTotal');

  
  
  