function Signal_to_BG_forMCPdosage
%% This script is to compare the signal to noise ratio for different MCP-mNG dosage.
% Scheme
% 1) Use LoadMS2Sets to load datasets with different MCP dosage, but taken
% with exactly the same imaging condition.
% which is : 48nm/pixel, 400nsec/pixel dwell time, 6 line accumulations

% 2) Calculate the signal to noise ratio
% How? : Use the first few frames for signal, then average to get averaged
% signal, and Offset averaged over all spots in the first few frames?

% 3) Calculate the MCP dosage
% How? : Calculate the averaged nuclear fluo -> best
% Alternative : Calculate the averaged offset (MS2 spots)

% 3) generate tuples of MCP-mNG dosage and S/N ratio

%% Load datasets
% Think about tweaking LoadMS2Sets to find files from designated folders.
filePath = 'E:\YangJoon\LivemRNA\Data\Dropbox\WormsMS2Results'

Prefix1 = '2020-01-22-sex1-50nm_6lineAcc_18sec_Series009'; % 2xMCP-mNG
Prefix2 = '2020-02-27-sex1-MS2-4xMCP-2'; %4xMCP-mNG
Prefix3 = '2020-03-10-sex1_LiAcc6_1'; 
Prefix4 = '2020-03-10-sex1_LiAcc6_2'; 
Prefix5 = '2020-03-10-sex1_LiAcc6_3'; 
%Prefix6 = '2020-03-10-sex1_LiAcc6_4'; 
Prefixes = {Prefix1, Prefix2, Prefix3, Prefix4, Prefix5};%, Prefix6};

%% for each dataset (Prefix)
% Use for loop to calculate the 
for j=1:length(Prefixes)
    % Load dataset
    Prefix = Prefixes{j};
    clear CompiledParticles
    load([filePath,filesep,Prefix,filesep,'CompiledParticles.mat'])
    load([filePath,filesep,Prefix,filesep,'Spots.mat'])
    % Calculate the average of spot fluorescence for the first few frames

    intArea = double(Spots(1).Fits(1).intArea);
    % Note that the offset should be multiplied with intArea (area which we
    % integrated the spot intensity)

    clear Fluo
    clear Fluo_mean
    clear FluoError
    clear Offset
    clear SNratio
    % 1. filter for particles existed for more than certain number of frames
    % (This is to filter ouf some false-positives)
    % 2. filter for the first few frames 
    frame_thresh = 2;
    frame_window = 1:frame_thresh;
    % 3. no negative fluo value

    k=1; % counter (particle)
    % Go through all particles to sort out the ones that can pass through the
    % filters described above.
    for i=1:length(CompiledParticles{1,1})
        if (length(CompiledParticles{1,1}(i).Frame) > frame_thresh) &&...
            ((CompiledParticles{1,1}(i).Frame(1) == 1)||...
            (CompiledParticles{1,1}(i).Frame(1) == 2)||...
            (CompiledParticles{1,1}(i).Frame(1) == 3))&&...
            sum(CompiledParticles{1,1}(i).Fluo(1:end)<0)==0 % no negative fluo at 1st/2nd frame

            Fluo = CompiledParticles{1,1}(i).Fluo;
            FluoError(k) = CompiledParticles{1,1}(i).FluoError; %Fluctuation of Offset (signal-to-noise)
            Fluo_mean(k) = nanmean(Fluo(frame_window)); % average spot fluo for the first 3 frames
        
            % offset for the nuc.fluo
            OffSet = CompiledParticles{1,1}(i).Off;
            Offset(k) = nanmean(OffSet(frame_window));
        
            % Signal-to-noise ratio
            SNratio(k) = Fluo_mean(k) / FluoError(k);
        
            % Signal-to-Background ratio
            SBratio(k) = Fluo_mean(k) / (Offset(k)*intArea);
            k=k+1;
        else
        end

    end
    
    % Some additional filter to consider outliers
    % 25-75% quantile?
    %percTile_25 = prctile(Fluo_mean,25);
    %percTile_75 = prctile(Fluo_mean,75);
    %Filter = (Fluo_mean > percTile_25) .* (Fluo_mean < percTile_75);
    % After plotting single spot fluorescence, I noticed that over 5500 is
    % somewhat outlier, most of which are false-positives such as membrane,
    % cytoplasmic junks, etc. Thus, I'll get rid of them.
    Filter = (Fluo_mean < 5000);
    Fluo_mean_filtered = Fluo_mean .*Filter;
    Fluo_mean_filtered(Fluo_mean_filtered==0) = []; % filter out 0s
    
    SpotFluo_filtered(j) = nanmean(Fluo_mean_filtered);
    SpotFluo_SEM_filtered(j) = nanstd(Fluo_mean_filtered,0)./sqrt(length(Fluo_mean_filtered));
    
    % Calculate the S/N ratio, SpotFluo, MCP-dosage(Offset)
    SNratio_mean(j) = nanmean(SNratio);
    SNratio_std(j) = nanstd(SNratio,0);
    SNratio_SEM(j) = SNratio_std(j) ./sqrt(k);

    SpotFluo(j) = nanmean(Fluo_mean);
    SpotFluo_SEM(j) = nanstd(Fluo_mean,0)./sqrt(length(Fluo_mean));
    
    % Individual spot fluo (averaged over 2-3 frames) saved to see how
    % variable they are
    ind_spotfluo{j} = Fluo_mean;
    ind_spotfluo_filtered{j} = Fluo_mean_filtered;
    
    SBratio_filtered = SBratio.*Filter;
    SBratio_filtered(SBratio_filtered ==0) = [];
    SBratio_filtered(j) = nanmean(SBratio_filtered);

    MCPdosage(j) = intArea * nanmean(Offset);
    MCPdosage_std(j) = intArea*nanstd(Offset,0);
end


%% Check individual spot fluo for different embryos 
% checking the variability (or possible mistakes in the quantification)
hold on
for i=1:length(ind_spotfluo)
    clear spotfluo
    clear embryo
    spotfluo = ind_spotfluo{i};
    MCP = MCPdosage(i)*ones(size(spotfluo));
    plot(MCP, spotfluo,'o')
end
xlim([0 40000])
xticks([0 10000 20000 30000 40000])
xlabel('MCP dosage')
ylabel('individual spot fluorescence (AU)')
StandardFigure(gcf,gca)

figPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Worms_MS2';
saveas(gcf,[figPath,filesep,'Signal_to_MCP_dosage_individualSpots.tif'])
saveas(gcf,[figPath,filesep,'Signal_to_MCP_dosage_individualSpots.pdf'])

%% Check individual spot fluo for different embryos 
% - after filtering with the maximum of 6000
% checking the variability (or possible mistakes in the quantification)
hold on
for i=1:length(ind_spotfluo_filtered)
    clear spotfluo
    clear embryo
    spotFluo_filtered = ind_spotfluo_filtered{i};
    MCP = MCPdosage(i)*ones(size(spotFluo_filtered));
    plot(MCP, spotFluo_filtered,'o')
end
xlim([0 40000])
ylim([0 5000])
xticks([0 10000 20000 30000 40000])
yticks([0 1000 2000 3000 4000 5000])
xlabel('MCP dosage(AU)')
ylabel('individual spot fluorescence (AU)')
StandardFigure(gcf,gca)

figPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Worms_MS2';
saveas(gcf,[figPath,filesep,'Signal_to_MCP_dosage_individualSpots_filtered.tif'])
saveas(gcf,[figPath,filesep,'Signal_to_MCP_dosage_individualSpots_filtered.pdf'])

%% Signal-to-Background
% Plot 
errorbar(MCPdosage, SpotFluo, SpotFluo_SEM,'o')

 xlim([0 40000])
% ylim([0 4])
 xticks([0 10000 20000 30000 40000])
xlabel('MCP dosage(AU)')
ylabel('Signal (AU)')
title('MS2 Signal for MCP dosage')
StandardFigure(gcf,gca)

figPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Worms_MS2';
% saveas(gcf,[figPath,filesep,'Signal_to_MCP_dosage.tif'])
% saveas(gcf,[figPath,filesep,'Signal_to_MCP_dosage.pdf'])

%% Signal-to-Background (after Quantile filter)
% Plot 
errorbar(MCPdosage, SpotFluo_filtered, SpotFluo_SEM_filtered,'o')

xlim([0 40000])
ylim([0 5000])
xticks([0 10000 20000 30000 40000])
yticks([0 1000 2000 3000 4000 5000])
xlabel('MCP dosage(AU)')
ylabel('Signal (AU)')
title('MS2 Signal for MCP dosage')
StandardFigure(gcf,gca)

figPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Worms_MS2';
saveas(gcf,[figPath,filesep,'Signal_to_MCP_dosage_filtered.tif'])
saveas(gcf,[figPath,filesep,'Signal_to_MCP_dosage_filtered.pdf'])


%% Signal-to-noise
errorbar(MCPdosage, SNratio_mean, SNratio_SEM,'o')

xlim([0 40000])
ylim([0 5])
xticks([0 10000 20000 30000 40000])
xlabel('MCP dosage(AU)')
ylabel('S/N ratio')
title('Signal-to-noise ratio for MCP dosage')
StandardFigure(gcf,gca)

figPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Worms_MS2';
saveas(gcf,[figPath,filesep,'SNratio_MCP_dosage.tif'])
saveas(gcf,[figPath,filesep,'SNratio_MCP_dosage.pdf'])

end