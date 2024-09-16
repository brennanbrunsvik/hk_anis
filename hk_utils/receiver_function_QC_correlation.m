function [rf, ac_keep]=receiver_function_QC_correlation(rf, tt, options)
    arguments
        rf
        tt
        options.acparms = [0.85, 1; 0.45, 40]; % Receiver function correlation parameters. [0.85, 1; 0.45, 40] means receiver function should have correlation of 0.85 in first 1 second, 0.45 in first 40 seconds. 
        options.plot_name = []; 
    end
    acparms = options.acparms; 

    % Remove those that have small zero-lag cross-correlation to other receiver functions. 
    % cross correlation parameters. These don't need to be perfect. It's just removing really bad receiver functions. 
    % First row: Handle bad parent pulse. Second is for whole receiver function. 
    % First collumn: Cutoff for noramlized average cross correlation to other
    % receiver functions. Second collumn: Time window (absolute value)
    
    fprintf("\nCleaning receiver functions with correlation of\n" + ...
        "%1.2f from -%1.2f to -%1.2f seconds and \n%1.2f from -%1.2f to %1.2f seconds\n", ...
        acparms(1,1), acparms(1,2), acparms(1,2), acparms(2,1), acparms(2,2), acparms(2,2))
    
    fprintf('\nStarting with %1.0f receiver functions.\n', size(rf,2))
    
    figure(61); clf; hold on; 
    subplot(2,1,1); hold on; 
    title('All')
    plot(tt, rf); 
    xlim([-5, 40]); 
    
    for iacremove = 1:2; % Loop over checking the first few seconds of the receiver function, and the whole receiver function. 
        nrf = size(rf,2); % Number of receiver functions. Changes with each iacremove
        ac = nan(nrf, nrf); % correlation matrix. 
        twin = abs(tt) < acparms(iacremove,2); % Window within which is the parent pulse, to cross-correlate.  
        ac_cutoff = acparms(iacremove,1); % Get parameters for this iacremove. 
        for irf = 1:nrf; % Loop over each receiver function
            for jrf = 1:nrf; % Again loop over each receiver function. This is fast enough and done not often so don't worry about taking advantage of symmetry
                rfi = rf(twin,irf); 
                rfj = rf(twin,jrf); 
                ac(irf, jrf) = rfi' * rfj ./ ...
                    ( sqrt(rfi' * rfi)*sqrt(rfj' * rfj) ); % cross-correlation. 
            end
        end
        acs = sum(ac) ./ length(ac); % Cross-correlation sum. Normalize it. 
        ac_keep = acs > ac_cutoff; % Only keep if there is high-enough average cross-correlation. 
        
        rf       = rf  (:,ac_keep     ); % Remove bad receiver functions. 
    end
    
    n_rf = size(rf,2)
    fprintf('\nEnding with %1.0f receiver functions.\n', n_rf)
    
    subplot(2,1,2); hold on; 
    title('Kept'); 
    plot(tt, rf); 
    xlim([-5, 40]); 

    if ~isempty(options.plot_name); 
        exportgraphics(gcf, options.plot_name); 
    end 

end
