function [finger1, finger2] = get_finger_pair(finger_pair_id)
    % Purpose %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get finger pair out of Thumb, Index, and Middle, using finger_pair_id  
    % Just for experiment 3 (Sweep experiment)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Input 
    %    finger_pair_id         : a -> (finger1, finger2) = (Thumb, Index)
    %                             b -> (finger1, finger2) = (Index, Middle)
    %                             c -> (finger1, finger2) = (Middle, Thumb)
    % Output 
    %    finger1                : Being applied frequency1  
    %    finger2                : Being applied frequency2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    switch finger_pair_id
        case 'a'
            finger1 = 'Thumb';
            finger2 = 'Index';
        case 'b'
            finger1 = 'Index';
            finger2 = 'Middle';
        case 'c' 
            finger1 = 'Middle';
            finger2 = 'Thumb';
        otherwise 
            disp('Set a correct id for getting the finger pair!');
    end

end
