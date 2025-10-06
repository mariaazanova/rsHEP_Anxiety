%% create processing directories
%Input: cell array of folder names

function [] = create_dirs(folderList, subjid)

if exist('subjid', 'var')
    for i = 1:length(folderList)
        if not(isfolder([folderList{i}, subjid]))
            mkdir([folderList{i}, subjid])
        end
    end
else
    
    for i = 1:length(folderList)
        if not(isfolder([folderList{i}]))
            mkdir([folderList{i}])
        end
    end
end

end