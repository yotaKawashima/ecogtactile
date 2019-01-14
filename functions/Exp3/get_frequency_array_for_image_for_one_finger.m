function fundamental_IM_array_for_image = get_frequency_array_for_image_for_one_finger(target_array_for_correlation, fundamental_array_for_correlation, IM_array_for_correlation)
    % Purpose %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Make frequency array for correspondence of group id to frequency
    % Just for experiment 3 (Sweep experiment)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Input 
    %    target_array_for_correlation      : Fundamental frequency array, whose frequency is allocated to the target finger    
    %    fundamental_array_for_correlation : Fundamental frequency array for correlation 
    %    IM_array_for_correlation          : IM frequency for correlation
    % Output 
    %    fundamental_IM_frequency_array_for_image   : Frequency array for image
    %    ----------
    %    dim = (faxis x condition x RGB)
    %    fundamental-> black, IM-> gray, background-> white
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Initialisation
    fundamental_IM_array_for_imageR = ones(size(fundamental_array_for_correlation,1), size(fundamental_array_for_correlation,2));
    fundamental_IM_array_for_imageG = ones(size(fundamental_array_for_correlation,1), size(fundamental_array_for_correlation,2));
    fundamental_IM_array_for_imageB = ones(size(fundamental_array_for_correlation,1), size(fundamental_array_for_correlation,2));
    
    % Set gray for fundamental frequency
    fundamental_IM_array_for_imageR(logical(fundamental_array_for_correlation)) = 180/256;
    fundamental_IM_array_for_imageG(logical(fundamental_array_for_correlation)) = 180/256;
    fundamental_IM_array_for_imageB(logical(fundamental_array_for_correlation)) = 180/256;
    
    % Set black for target fundamental frequency
    fundamental_IM_array_for_imageR(logical(target_array_for_correlation)) = 0;
    fundamental_IM_array_for_imageG(logical(target_array_for_correlation)) = 0;
    fundamental_IM_array_for_imageB(logical(target_array_for_correlation)) = 0;
    
    % Set pink for IM frequency
    if nargin > 3
        fundamental_IM_array_for_imageR(logical(IM_array_for_correlation)) = 240/256;
        fundamental_IM_array_for_imageG(logical(IM_array_for_correlation)) = 0/256;
        fundamental_IM_array_for_imageB(logical(IM_array_for_correlation)) = 120/256;
    end
    
    fundamental_IM_array_for_image = cat(3, fundamental_IM_array_for_imageR, fundamental_IM_array_for_imageG, fundamental_IM_array_for_imageB); 
end
