function [ascending_frequency, descending_frequency] = make_frequency_array()
    % Purpose %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Make frequency array for correspondence of group id to frequency
    % Just for experiment 3 (Sweep experiment)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Input 
    %    None
    % Output 
    %    ascending_frequency    : Frequency array for ascending condition 
    %    descending_frequency   : Frequency array for descending condition
    %    ----------
    %    dim = (21 x 2)
    %    [[frequency1, frequency2],
    %     [frequency1, frequency2],
    %              ...
    %     [frequency1, frequency2]]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Ascending
    ascending_frequency = zeros(21,2);
    frequency1 = 8 ;
    frequency2 = 10 ;
    for f_id = 1:21
        % Store frequency
        ascending_frequency(f_id,:) = [frequency1, frequency2];

        % Update frequency
        frequency1 = frequency1 + 1;
        frequency2 = frequency2 + 2;
    end % f_id

    % Descending 
    descending_frequency = zeros(21,2);
    frequency1 = 50 ;
    frequency2 = 48 ;
    for f_id = 1:21
        % Store frequency
        descending_frequency(f_id,:) = [frequency1, frequency2];

        % Update frequency
        frequency1 = frequency1 - 1;
        frequency2 = frequency2 - 2;
    end % f_id
end
