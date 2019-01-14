function frequency_array_for_correlation = get_frequency_array_for_correlation(faxis, frequency_array)
    % Purpose %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Make frequency array for correspondence of group id to frequency
    % Just for experiment 3 (Sweep experiment)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Input 
    %    faxis                  : Frequency array from FFT
    %    frequency_array        : Frequency array for trial condition
    % Output 
    %    frequency_array_for_correlation   : Frequency array for plot
    %    ----------
    %    dim = (faxis x condition)
    %    if frequency on, then assign 1 to the element of the array
    %    else (frequency off), then assign 0 to the element of the array
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Initialisation 
    frequency_array_for_correlation = zeros(length(faxis), 21);
    
    for condition_id = 1:21
        for frequency_id = 1:size(frequency_array, 2)
            % Get frequency
            freq_for_one_finger = frequency_array(condition_id, frequency_id);
            
            % Get index
            [~,idxF]=findclosest(faxis,freq_for_one_finger);
            
            % Assign 1 to the element
            frequency_array_for_correlation(idxF, condition_id) = 1;
            
        end
    end
    
end
