function checkBleaching_wholeImage(Prefix)
% This script is for checking the degree of bleaching using the whole
% fluorescence in the field of view.
% For example, let's say there's a movie with xyzt.
% I can use the mean (median) of fluorescence of the whole field of view
% over z-stacks/x-y plane at each time point, then plot it to check the
% bleaching.
% If the fluorescence doesn't change that much, then it means that the
% bleaching is negligible.

% To implement
% 1) I need to save the results with the specific imaging
% condition, which I can pull from the metadata.
% 2) 

% input : Prefix
% output : fluo (with SD) vs time
%% Load the dataset using Prefix

% Define the Prefix (either here or as an input argument)
Prefix = '2020-02-27-sex1-MS2-4xMCP-1';

% Read images from the Pre-Processed folder
filePath = 'E:\Worms\Data\PreProcessedData';

% Count the # of frames using the number of His** files.
D = [filePath, filesep, Prefix];
List = dir(D);
filenames = {List.name};
numFrames = sum( ~cellfun(@isempty, strfind(filenames, 'His')) )

for i=1:numFrames
    if i<10
        index = strcat('_00',num2str(i));
    else
        index = strcat('_0',num2str(i));
    end
    
    % Second, read all the z-slices
    numZslices =  sum( ~cellfun(@isempty, strfind(filenames, index))) -1; % -1 is for His image.
    for j=2:numZslices-1
        % Go over specific frame, and z-slice one by one
        tIndex = (~cellfun(@isempty, strfind(filenames, index)));
        if j<10
            zindex = strcat('z0',num2str(j));
        else
            zindex = strcat('z',num2str(j));
        end
        
        zIndex = (~cellfun(@isempty, strfind(filenames, zindex)));
        
        Index = find(tIndex.*zIndex);
        
        Image(:,:,j) = imread([filePath, filesep, Prefix, filesep, List(Index).name]);
        Intensity_sum(j) =sum(sum(Image(:,:,j)));
        
    end
    %Intensity_median(i) = nanmedian(Image,'all'); % median for all x,y,z direction
    %Intensity_mean(i) = nanmean(Image,'all'); % mean for all x,y,z direction
    %Intensity_std(i) = nanstd(double(Image),0,'all');
    
    Intensity_mean(i) = nanmean(Intensity_sum);
    Intensity_SEM(i) = nanstd(Intensity_sum,0)./sqrt(numZslices-2);
    %Image_Zprojected(:,:,i) = nanmax(Image,[],3);
  
    
end

% Note that we need to define the channel

%% Plot the mean intensity over x,y,z image at each time point
% Load compiledparticles
% filePath = 'E:\YangJoon\LivemRNA\Data\Dropbox_old\WormsMS2Results'
% load([filePath,filesep,Prefix,filesep,'CompiledParticles.mat'])

errorbar(ElapsedTime, Intensity_mean / max(Intensity_mean), Intensity_SEM / max(Intensity_mean))
ylim([0 1.2])


%% 
hold on

errorbar(Time_LineAcc1, Intensity_mean_LinAcc1 / max(Intensity_mean_LinAcc1),...
            Intensity_SEM_LinAcc1 / max(Intensity_mean_LinAcc1))
errorbar(Time_LineAcc6, Intensity_mean_LinAcc6 / max(Intensity_mean_LinAcc6),...
            Intensity_SEM_LinAcc6 / max(Intensity_mean_LinAcc6))
ylim([0 1.2])

legend('Line Accumulation 1','Line Accumulation 6')
title('Bleaching check')
xlabel('time (min)')
ylabel('normalized pixel intensity(AU)')

StandardFigure(gcf,gca)

% save the figure
figPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Worms_MS2';
saveas(gcf,[figPath,filesep,'BleachingTest_48nm_pixel_400nsec_LineAccumulations.tif'])
saveas(gcf,[figPath,filesep,'BleachingTest_48nm_pixel_400nsec_LineAccumulations.pdf'])

%% Create a structure to compile the results per embryo (condition/dosage)
Bleaching_check{1,1} = 'Prefix';
Bleaching_check{1,2} = 'Intensity_mean (x-y) over z, @ t';
Bleaching_check{1,3} = 'Intensity_SEM';
Bleaching_check{1,4} = 'ElpasedTime';

Bleaching_check{2,1} = '2020-02-27-sex1-MS2-4xMCP-1';
Bleaching_check{2,2} = Intensity_mean_LinAcc1;
Bleaching_check{2,3} = Intensity_SEM_LinAcc1;
Bleaching_check{2,4} = Time_LineAcc1;

Bleaching_check{3,1} = '2020-02-27-sex1-MS2-4xMCP-2';
Bleaching_check{3,2} = Intensity_mean_LinAcc6;
Bleaching_check{3,3} = Intensity_SEM_LinAcc6;
Bleaching_check{3,4} = Time_LineAcc6;


%% Save the processed data

end