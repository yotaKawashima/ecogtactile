%% Plot intermodulation

clear all;
close all;

localdef;

f1 = 50;
f2 = 48;

order_struct = struct();

nOrder = 8; 

order_struct = struct();
order_struct.order = 1:1:nOrder;
start_id = zeros(1, nOrder);

counter_set = 0;
    
for i_order = 1:nOrder
    if i_order == 1
        coefficience_set = [0 1];
        coefficience_set = cat(1, coefficience_set, [1 0]);
        start_id(i_order) = 1;
        counter_set = 2;
    else
        for element_1 = 1:i_order-1
            element_2 = i_order - element_1;
            coefficience_set = cat(1, coefficience_set, [element_1 element_2]);
            counter_set = counter_set + 1;
            if element_1 == 1
                start_id(i_order) = counter_set;
            end 
        end % element_1
    end % if i_order
end % i_order

order_struct.set_start_id = start_id;
order_struct.sets = coefficience_set;

%% Plot ims

fig_1 = figure();

for i_order  = 1:nOrder
    
    y_value = (11 - i_order) ;
    
    start_id = order_struct.set_start_id(i_order);
    
    if i_order ~= nOrder
        end_id = order_struct.set_start_id(i_order + 1) -1; 
    else
        end_id = size(order_struct.sets, 1);
    end 
    
    ims = zeros(size(start_id:end_id,1)*2, 1);
    
    counter_im = 1;
    
    if i_order == 1
        for i_set = start_id:end_id 
            set_now = order_struct.sets(i_set,:);
            ims(counter_im) = f1*set_now(1) + f2*set_now(2);
            
            figure(fig_1);
            hold on;
            line([ims(counter_im) ims(counter_im)], [0 y_value], 'Color', 'k');
            counter_im = counter_im + 1;
        end % i_set
    else
        for i_set = start_id:end_id 
            set_now = order_struct.sets(i_set,:);
            for sign = 1:2
                if sign == 1
                    ims(counter_im) = f1*set_now(1) + f2*set_now(2);
                elseif sign == 2
                    ims(counter_im) = abs(f1*set_now(1) - f2*set_now(2));
                end

                figure(fig_1);
                hold on;
                line([ims(counter_im) ims(counter_im)], [0 y_value], 'Color', 'k');
                counter_im = counter_im + 1;

            end % sign
        end % i_set
    end % i_order
end % i_order

grid on;
ylim([0,10]);
yticklabels(["","","","8","7","6","5","4", "3", "2", "1"]);
title('Intermodulation spectrum');
xlabel('Frequency [Hz]');
ylabel('nth-order IM [-]');
format_fig;