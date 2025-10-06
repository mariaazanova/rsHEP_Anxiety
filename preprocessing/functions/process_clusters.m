%% Define helper function
%% plot significant clusters
function process_clusters(stat, clustfield, labelfield, sign_prefix, dataFT, hp_settings, baseline, ic_cutoffs, hp, bl, ic, timecourse_dir, topo_dir)
    if isfield(stat, clustfield) && ~isempty(stat.(clustfield))
        pvals = [stat.(clustfield).prob];
        sigIdx = find(pvals < 0.05);
        
        for iClust = sigIdx
            clustMask = (stat.(labelfield) == iClust);
            
            % Find channels that are part of the cluster at any time point
            sigChans = any(clustMask, 2);  % logical index of significant channels
            
            % Average across significant channels for all time points
            clusterWave = squeeze(mean(dataFT.individual(:,sigChans,:), [1 2]));

            % Plot and save time course
            h1 = figure; hold on;
            
            % Find significant time window
            sigTime = any(clustMask, 1);
            timeIdx = find(sigTime);
            
            % Add gray patch for significant time window
            yLims = [min(clusterWave) max(clusterWave)];
            patch([dataFT.time(timeIdx(1)) dataFT.time(timeIdx(end)) ...
                   dataFT.time(timeIdx(end)) dataFT.time(timeIdx(1))], ...
                  [yLims(1) yLims(1) yLims(2) yLims(2)], ...
                  [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
            
            % Plot time course on top
            plot(dataFT.time, clusterWave, 'LineWidth', 2);
            xlabel('Time (s)'); ylabel('Amplitude (\muV)');
            title(sprintf('%s Cluster %d Average Timecourse (p=%0.3f)', ...
                sign_prefix, iClust, stat.(clustfield)(iClust).prob));
            axis tight; grid on;
            
            % Save time course plot
            fname_time = sprintf('timecourse_hp_%0.1f_bl_%s_ic_%0.1f_%s_clust%d.png', ...
                hp_settings(hp), mat2str(baseline{bl}), ic_cutoffs(ic), ...
                lower(sign_prefix), iClust);
            saveas(h1, fullfile(timecourse_dir, fname_time));
            close(h1);
            
            % Create topography by averaging channels over significant time window
            timeMask = any(clustMask, 1);
            meanVals = squeeze(mean(mean(dataFT.individual(:,:,timeMask), 1), 3)); 

            topoData = [];
            topoData.avg = meanVals;
            topoData.label = dataFT.label;
            topoData.time = 0;        % dummy single time point
            topoData.dimord = 'chan_time';

            h2 = figure;
            cfgp = [];
            cfgp.layout = 'easycapM11.lay';
            cfgp.xlim = [0 0];
            ft_topoplotER(cfgp, topoData);
            title(sprintf('%s Cluster %d Topography (p=%0.3f)', ...
                sign_prefix, iClust, stat.(clustfield)(iClust).prob));

            % Save topography plot
            fname_topo = sprintf('topo_hp_%0.1f_bl_%s_ic_%0.1f_%s_clust%d.png', ...
                hp_settings(hp), mat2str(baseline{bl}), ic_cutoffs(ic), ...
                lower(sign_prefix), iClust);
            saveas(h2, fullfile(topo_dir, fname_topo));
            close(h2);
        end
    end
end
