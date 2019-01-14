function fundamental_IM_array_for_image = get_frequency_array_for_image_for_one_finger_4IM(target_array_for_correlation, fundamental_array_for_correlation, IM1_array_for_correlation, IM2_array_for_correlation, IM3_array_for_correlation, IM4_array_for_correlation)
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
    
    % Set pink for IM frequency
    fundamental_IM_array_for_imageR(logical(IM1_array_for_correlation)) = 137/256;
    fundamental_IM_array_for_imageG(logical(IM1_array_for_correlation)) = 15/256;
    fundamental_IM_array_for_imageB(logical(IM1_array_for_correlation)) = 80/256;
    
    fundamental_IM_array_for_imageR(logical(IM2_array_for_correlation)) = 256/256;
    fundamental_IM_array_for_imageG(logical(IM2_array_for_correlation)) = 24/256;
    fundamental_IM_array_for_imageB(logical(IM2_array_for_correlation)) = 69/256;
    
    fundamental_IM_array_for_imageR(logical(IM3_array_for_correlation)) = 245/256;
    fundamental_IM_array_for_imageG(logical(IM3_array_for_correlation)) = 144/256;
    fundamental_IM_array_for_imageB(logical(IM3_array_for_correlation)) = 178/256;
    
    fundamental_IM_array_for_imageR(logical(IM4_array_for_correlation)) = 253/256;
    fundamental_IM_array_for_imageG(logical(IM4_array_for_correlation)) = 229/256;
    fundamental_IM_array_for_imageB(logical(IM4_array_for_correlation)) = 237/256;
        
    % Set black for target fundamental frequency
    fundamental_IM_array_for_imageR(logical(target_array_for_correlation)) = 0;
    fundamental_IM_array_for_imageG(logical(target_array_for_correlation)) = 0;
    fundamental_IM_array_for_imageB(logical(target_array_for_correlation)) = 0;
    
    fundamental_IM_array_for_image = cat(3, fundamental_IM_array_for_imageR, fundamental_IM_array_for_imageG, fundamental_IM_array_for_imageB); 
end
