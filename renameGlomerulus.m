% Rename a glomerulus in a specified folder
clear all
close all

folderToRename='';
labelName='manualLabels.mat';

newGlomName='DC2';
oldGlomName='DC1';

cd(folderToRename)

currfolders=dir(pwd);
currfolders=currfolders(3:end);

for i=1:length(currfolders)
    cd(currfolders(i).name)
    
    currfiles=dir(pwd);
    for j=1:length(currfiles)
        if strfind(currfiles(j).name,labelName)
            currName=currfiles(j).name;
            break
        end
    end
    
    load(currName)
    
    yesrename=0;
    % first check to see if oldGlomName exists. if not, do not try to
    % rename
    for j=1:length(clusterLabels)
        if strcmp(clusterLabels{j},oldGlomName)
            yesrename=1;
        end
    end
    
    if yesrename
        % first remove any labels referencing the new glom name
        for j=1:length(clusterLabels)
            if strcmp(clusterLabels{j},newGlomName)
                clusterLabels{j}=num2str(j);
            end
        end
        
        % rename oldGlomName to newGlomName
        for j=1:length(clusterLabels)
            if strcmp(clusterLabels{j},oldGlomName)
                clusterLabels{j}=newGlomName;
            end
        end
        
        save(currName,'clusterLabels','savename')
    end
    cd ..
end