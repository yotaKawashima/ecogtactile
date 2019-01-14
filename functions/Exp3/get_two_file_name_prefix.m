function [file_name_prefix1, file_name_prefix2] = get_two_file_name_prefix(file_name_prefixs, sweep_id, frequency_id, finger)
    % Purpose %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get two file name prefix for the same finger, the same frequency, and the same sweep type   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Input 
    %    file_name_prefixs          : all file name prefixs (cell array) 
    %    sweep_id                   : Sweep_id (Ascending = 1, Descending = 2)
    %    frequency_id               : Frequency_id (frequency1 = 1, frequency2 = 2)
    %    finger                     : Finger(Thumb/Index/Middle)
    % Output 
    %    file_name_prefix1/2        : File name prefix 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    switch frequency_id 
        case 1
            switch finger
                case 'Thumb'
                    file_name_prefix1 = file_name_prefixs{(6*(sweep_id -1) + 1),1}; % (Thumb, Index)  = (frequency1, frequency2)
                    file_name_prefix2 = file_name_prefixs{(6*(sweep_id -1) + 6),1}; % (Thumb, Middle) = (frequency1, frequency2)
                case 'Index'
                    file_name_prefix1 = file_name_prefixs{(6*(sweep_id -1) + 2),1}; % (Index, Middle) = (frequency1, frequency2)
                    file_name_prefix2 = file_name_prefixs{(6*(sweep_id -1) + 4),1}; % (Index, Thumb)  = (frequency1, frequency2)
                case 'Middle'
                    file_name_prefix1 = file_name_prefixs{(6*(sweep_id -1) + 3),1}; % (Thumb, Index)  = (frequency1, frequency2)
                    file_name_prefix2 = file_name_prefixs{(6*(sweep_id -1) + 5),1}; % (Thumb, Middle) = (frequency1, frequency2)
            end

        case 2
            switch finger
                case 'Thumb'
                    file_name_prefix1 = file_name_prefixs{(6*(sweep_id -1) + 3),1}; % (Middle, Thumb) = (frequency1, frequency2)
                    file_name_prefix2 = file_name_prefixs{(6*(sweep_id -1) + 4),1}; % (Index,  Thumb) = (frequency1, frequency2)
                case 'Index'
                    file_name_prefix1 = file_name_prefixs{(6*(sweep_id -1) + 1),1}; % (Thumb,  Index) = (frequency1, frequency2)
                    file_name_prefix2 = file_name_prefixs{(6*(sweep_id -1) + 5),1}; % (Middle, Index) = (frequency1, frequency2)
                case 'Middle'
                    file_name_prefix1 = file_name_prefixs{(6*(sweep_id -1) + 2),1}; % (Index, Middle) = (frequency1, frequency2)
                    file_name_prefix2 = file_name_prefixs{(6*(sweep_id -1) + 6),1}; % (Thumb, Middle) = (frequency1, frequency2)
            end
    end

end
