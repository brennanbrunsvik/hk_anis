function [rf, ac_keep]=receiver_function_QC_correlation(rf, tt, options)
    arguments
        rf
        tt
        options.acparms = [0.85, 1; 0.45, 40]; % Receiver function correlation parameters. [0.85, 1; 0.45, 40] means receiver function should have correlation of 0.85 in first 1 second, 0.45 in first 40 seconds. 
        options.actime_cluster = 40; % If using clustering, then do so for receiver functions over this time span. 
        options.plot_name = []; 
        options.noise_removal_style = 'cluster'; % "cluster" or "average_dist". clustering (with dbscan), or average dissimilarity to other receiver functions. 
    end
    acparms = options.acparms; 
    % rf: Receiver functions. ntt x nrf. 
    % tt: Time. ntt x 1. 
    % Remove those that have small zero-lag cross-correlation to other receiver functions. 
    % cross correlation parameters. These don't need to be perfect. It's just removing really bad receiver functions. 
    % If using "average_dist": 
    %     First row: Handle bad parent pulse. Second is for whole receiver function. 
    %     First collumn: Cutoff for noramlized average cross correlation to other
    %     receiver functions. Second collumn: Time window (absolute value)
    % If using "cluster": 
    %     No input parameters needed as of 2024/09/16

    if size(tt,1) ~= size(rf,1); 
        error('brb2024/09/16: Input sizes are not correct. rf: %dx%d, tt: %dx%d', size(rf, 1), size(rf, 2), size(tt, 1), size(tt, 2) );
    end

    fprintf('\nStarting to clean receiver functions.\n')
    fprintf('\nStarting with %1.0f receiver functions.\n', size(rf,2))

    figure(61); clf; hold on; 
    subplot(2,1,1); hold on; 
    title('All')
    plot(tt, rf); 
    xlim([-5, 40]); 

    if strcmp(options.noise_removal_style, 'cluster'); 
        fprintf('\nCleaning receiver functions using clustering. \n')
        %% Clustering approach. 
   
        nrf = size(rf, 2); 
        twin = abs(tt) < options.actime_cluster;
    
        % Precompute the norms of each receiver function within the twin window
        rfnorm = sqrt(sum(rf(twin, :) .^ 2));
    
        % Vectorized cross-correlation computation using matrix operations
        ac = (rf(twin, :)' * rf(twin, :)) ./ (rfnorm' * rfnorm);
    
        % Convert from similarity to distance metric. 
        dist = 1 - ac; 
        for ii = 1:size(dist,1); 
            dist(ii,ii) = 0; 
        end 

    
        % Determine clustering parameters. 
        forstd = dist(~eye(size(dist))); % Don't include the arbitarary 0 in distance mean calculation
        epsilon = 2*std(forstd); % Standard deviation of similarities. If a receiver function cannot be reached within 1 std from the cluster, it seems reasonable to call it noise. 
        minpts = ceil(.1 * nrf); % 30% of number of receiver functions. From Petruska and Eilon 2022 GGG. We just get one cluster this way, which is good for separating nosie.  
        
        [dbres] = dbscan(dist, epsilon, minpts, 'Distance', 'precomputed'); % DBSCAN clustering. 
        unique_dbres = unique(dbres); % Number of clusters. -1 is noise. 
        ac_keep = dbres ~= -1; % Keep everything that was not marked as noise. 
        ratio_keep = sum(ac_keep) ./ length(ac_keep); 

        fprintf('Number of unique clusters, including noise: %1.0f.\n', length(unique_dbres))


        rf = rf(:, ac_keep);        

    elseif strcmp(options.noise_removal_style, 'average_dist'); 
        fprintf('\nCleaning receiver functions using simple distance metric. \n')
    
        fprintf("\nCleaning receiver functions with correlation of\n" + ...
            "%1.2f from -%1.2f to -%1.2f seconds and \n%1.2f from -%1.2f to %1.2f seconds\n", ...
            acparms(1,1), acparms(1,2), acparms(1,2), acparms(2,1), acparms(2,2), acparms(2,2))
    
        fprintf('\nStarting with %1.0f receiver functions.\n', size(rf,2))
    
        figure(61); clf; hold on; 
        subplot(2,1,1); hold on; 
        title('All')
        plot(tt, rf); 
        xlim([-5, 40]); 
    
        for iacremove = 1:2
            nrf = size(rf, 2); 
            ac = nan(nrf, nrf);
            twin = abs(tt) < acparms(iacremove, 2);
            ac_cutoff = acparms(iacremove, 1);
    
            % Precompute the norms of each receiver function within the twin window
            rfnorm = sqrt(sum(rf(twin, :) .^ 2));
    
            % Vectorized cross-correlation computation using matrix operations
            ac = (rf(twin, :)' * rf(twin, :)) ./ (rfnorm' * rfnorm);
    
            % Cross-correlation sum and normalization
            acs = sum(ac, 2) ./ nrf;
    
            % Remove bad receiver functions
            ac_keep = acs > ac_cutoff;
            rf = rf(:, ac_keep);
        end
    end

    n_rf = size(rf,2); 
    fprintf('\nEnding with %1.0f receiver functions.\n', n_rf)
    figure(61); 
    subplot(2,1,2); hold on; 
    title('Kept'); 
    plot(tt, rf); 
    xlim([-5, 40]); 

    if ~isempty(options.plot_name); 
        exportgraphics(gcf, options.plot_name); 
    end


end
