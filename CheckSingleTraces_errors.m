%% Checking single MS2 trace with error bars
function CheckSingleTraces_errors(varargin)
%% Load the datasets
% CompiledParticles
% Define the directory
% Dropbox folder etc.
% Prefix = varargin{1};
% Prefix = '2019-12-23-Sdc2-24xMS2-pMex5-MCP-GFP-200nmPixel'; % Particle23
% Prefix = '2019-12-20-sex1-24xMS2-pMex5-MCP-GFP-Series024_70nmPixel'; %Particle 20 and 25 are in the same nucleus
% Prefix = '2020-01-06-sex1-MS2-pMex5-MCP-GFP-LiAcc3-1'; % Particle #73, 74
% Prefix = '2020-01-22-sex1-50nm_6lineAcc_18sec_Series009'; % Particle 1&3, 5&7 (two nuclei)
%Prefix = '2020-02-27-sex1-MS2-4xMCP-1'; % Line Accumulation = 1, 400nsec, 48nm/pixel
% Prefix = '2020-02-27-sex1-MS2-4xMCP-2'; % Line Accumulation = 6, 400nsec, 48nm/pixel
Prefix = '2020-03-10-sex1_LiAcc6_1'
filePath = 'E:\YangJoon\LivemRNA\Data\Dropbox\WormsMS2Results';


% Prefix = '2019-08-14-hbP2-r2_close_MS2V5-lacZ-2';
% filePath = 'E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient';

load([filePath, filesep, Prefix, filesep, 'CompiledParticles.mat']);
load([filePath, filesep, Prefix, filesep, 'FrameInfo.mat']);
%load([filePath, filesep, Prefix, filesep, 'Spots.mat']);

CompiledParticles = CompiledParticles{1,1};
%% Pick a good particle that I manually curated
hold on

%Prefix = '2020-02-27-sex1-MS2-4xMCP-1'; % Line Accumulation = 1, 400nsec, 48nm/pixel
%[3,6,13,18]

particleIndex = 22; 
% particleIndex = 59;

Fluo = CompiledParticles(particleIndex).Fluo;
% Note that the error is the fluctuation of the Offset, 
% and it should be multiplied with the
% integration area (intArea) to match the dimension.\

% The FluloError field defined in GetParticleTrace.m is that value, already
% multiplied with the area of integration.

%spots_firstFrame = Spots(1).Fits;
%intArea = spots_firstFrame(1).intArea;
Offset = CompiledParticles(particleIndex).Off; %.* double(intArea);
% Fluo Error field is already multiplied with the integration Area, and
% also it's from the fluctuation of the offset. (Time-independent)
FluoError = CompiledParticles(particleIndex).FluoError;
FluoError = ones(size(Fluo)) * FluoError;
% time frames
Frames = CompiledParticles(particleIndex).Frame;

errorbar(ElapsedTime(Frames), Fluo, FluoError)
%plot(ElapsedTime(Frames), Offset)
title('Single MS2 trace')
xlabel('Time (min)')
ylabel('Fluorescence (AU)')

% Save 
end