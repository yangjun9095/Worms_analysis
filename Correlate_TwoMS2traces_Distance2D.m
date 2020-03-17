% Script for worms MS2 data
% Generating plots for two spots from the same nucleus

% In future, we need to pull the info about the nucleus where the MS2 spots (particles)
% are associated with. But, for now, I'll pick spots based on
% "CheckParticleTracking" indices.

%% Load the dataset (Particles, Spots, and FrameInfo)
filePath = 'E:\YangJoon\LivemRNA\Data\Dropbox\WormsMS2Results\2019-12-20-sex1-24xMS2-pMex5-MCP-GFP-Series024_70nmPixel'

load([filePath, filesep, 'FrameInfo.mat'])
load([filePath, filesep, 'Particles.mat'])
load([filePath, filesep, 'CompiledParticles.mat'])
load([filePath, filesep, 'Spots.mat'])

%% Particles : #20, and #25 in this case.
CompiledParticles = CompiledParticles{1,1};
frames1 = CompileParticles(20).Frames;
fluo1 = CompileParticles(20).Fluo;

frames2 = CompileParticles(25).Frames;
fluo2 = CompileParticles(25).Fluo;
%% 