%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Pipeline for dual EEG recorings, final study with drifting metronomes 
% Data collection started on 27/01/2020 in Ghent.
%
%
% Import in the same workspace a dual EEG recording; as configurations are
% the same for both subjects, run loops for subi = 1:2 in each section; use
% parfor only when computations are heavier (e.g. ICA), i.e. when it's worth 
% the overhead. When windows for visual inspection are to be opened, use
% for loops.
% Initialize cells for the dyad's dataset at each step before every for/parfor
% loop. Cells are better because sizes might differ and are necessary when
% parallel computing is used, since they might not store simultaneously
% resulting in different array sizes.
%
% Blocks need to be run one at a time: at the end of every step, a new
% dataset will be saved in the workspace. e.g.: contData = continuous
% data; filtData = filtered data; trialData = data segmented into trials.


% Mattia Rosso (Ghent, 11/01/2021)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
clc

% Set Fieldtrip path
%https://www.fieldtriptoolbox.org/faq/should_i_add_fieldtrip_with_all_subdirectories_to_my_matlab_path/
restoredefaultpath
addpath '/Users/mattiaipem/Documents/MATLAB/fieldtrip-20201205' % put your own fieldtrip path
ft_defaults


% Set your own path
path = '/Users/mattiaipem/Desktop/Hyperscanning_data/EEG/eegData/eeg_raw';
cd(path);


% Pick dyad to preprocess
pair = 1;

% Order of conditions, and design settings
condorder = [1 2 3 4 ; 4 3 1 2 ; 1 2 4 3 ; 1 2 3 4 ; 2 1 3 4 ; 4 3 1 2 ; 3 4 2 1 ; ...
             3 4 1 2 ; 4 3 1 2 ; 4 3 2 1 ; 2 1 4 3 ; 4 3 1 2 ; 2 1 4 3 ; 4 3 1 2];
blocks  =   [1:10 ; 11:20 ; 21:30 ; 31:40];
ntrials =   40;

% Channel division per subjects (total amount of channels = 128)
splitChan = {[1:64] , [65:128]}; %might experiment with PARFOR, for heavier steps (you can test it with copies of the same pair, simulating having a real sample)
% Remove lables to make channels recognizable by FieldTrip; re-label them after ICA artifact removal
prefix = {'P1_' , 'P2_'}; 

% Set sampling rate
srate = 1000;

% Find trigger's code and samples
events = ft_read_event(['Coordination_' num2str(pair) '.vhdr']);
% For mistakes with triggers have a look at ft_recodeevent(cfg, data). 

% Epoching in trials
% Refresh configuration at each step
cfg=[];
cfg.dataset = ['Coordination_' num2str(pair) '.vhdr']; %header to read
cfg.trialdef.eventtype = 'Stimulus'; %ask FieldTrip, setting '?'
if pair == 5 %Mistake: I forgot to change the trigger: condition 2 should was labelled as baseline (didn't change the trigger cable) 
    cfg.trialdef.eventvalue = {'s1' , 's2' , 's4' , 's8' , 's16'};
else
    cfg.trialdef.eventvalue = {'s1' , 's2' , 's4' , 's8'};
end
cfg.trialdef.prestim  = 0;   %I will pick the baseline excluding count-in time
cfg.trialdef.poststim = 39.008;   %in sec
cfg.trialdef.ntrials  = inf;   %all trials
cfg.trialfun = 'ft_trialfun_general';
% Define trials in dedicated configuration; needed before ft_preprocessing,
% but applied later
[cfg_trial] = ft_definetrial(cfg);

% Big cell arrays
contData  = {[]};
filtData  = {[]};
trialData = {[]};
eog = {[]}; % to store EOG channel

% Loop over subjects within the pair
for subi = 1:2 
    
    % Exceptionally, here we configure inside the loop as the subjects have
    % different subsets of channels
    cfg = [];
    cfg.continous = 'yes';
    cfg.dataset = ['Coordination_' num2str(pair) '.vhdr'];
    cfg.channel = splitChan{subi};
    cfg.demean = 'yes';   %it means indeed removing the DC offset (subtracting mean in the the frequency domain)
    cfg.implicitref = 'CPz';
    % recover CPz; comment for checking the position
    % of the electrode; reflect on how should I do with sub2_CPz
    
    % Process continuous data
    contData{subi} = ft_preprocessing(cfg);
    % remove prefix from labels
    contData{subi}.label = strrep(contData{subi}.label , prefix{subi} , '') ;
    
    % Note: IIR Butterworth filters are default; order is 6 by default (see low-level FT functions)
    cfg = [];
    cfg.bsfilter = 'yes'; %notch
    cfg.bsfreq = [49 51; 99 101; 149 151];
    cfg.hpfilter = 'yes'; %high-pass filter for slow drifts
    cfg.hpfreq = 1;
    cfg.lpfilter = 'yes'; %low-pass filter for hi-frex muscular activity
    cfg.lpfreq = 45;
    filtData{subi} = ft_preprocessing(cfg , contData{subi});
    
    % Re-arrange electrodes to match layout
    % Store EOG in separate cell; to further remove from channels
    eog{subi} = filtData{subi}.trial{1}(32,:);
    % exclude here EOG
    cfg              = [];
    cfg.channel      = {'all' '-EOG'};
    filtData{subi}   = ft_selectdata(cfg, filtData{subi});
    % move CPz in its position (49), and Oz in 31 
    tempLabel = [filtData{subi}.label(1:30);filtData{subi}.label(end-1);filtData{subi}.label(31:47);filtData{subi}.label(end);filtData{subi}.label(48:end-2)];
    tempData = [filtData{subi}.trial{1}(1:30,:);filtData{subi}.trial{1}(end-1,:);filtData{subi}.trial{1}(31:47,:);filtData{subi}.trial{1}(end,:);filtData{subi}.trial{1}(48:end-2,:)];
    filtData{subi}.label = tempLabel;
    filtData{subi}.trial{1} = tempData;
    clear tempLabel;
    clear tempData;
    
    % Segmenting function with configuration and continuous filtered data as input argument
    trialData{subi} = ft_redefinetrial(cfg_trial , filtData{subi});

%{
    % For inspection
    cfg.viewmode = 'vertical';
    ft_databrowser(cfg,trialData{subi});
    pause;
%}
end

% Correct the 90deg rotation, using filtData as imput (eog removed, CPz included);
% from now on, use 'layout'
cfg = [];
cfg.layout    = 'standard_waveguard64.elc';
cfg.rotate      = 90;
layout = ft_prepare_layout(cfg , filtData{1});

% %Uncomment for inspecting layout
% cfg = [];
% cfg.layout    = layout;
% ft_layoutplot(cfg, filtData{1});  %data before channels rejection
% 

  
%% Visual rejection trials and chans

%If you have a long recording, you will notice things tend to get worse
%towards the end, due to gel drying and electrodes losing contact
%NOTE: in this case, I have a lot of electrodes but few trials... don't
%sacrifice trials unless it's strictly necessary

%ALWAYS CHECK for electrodes with variance 0 which are not the reference!! (i.e., 'flat'). They MUST be
%removed, or the ICA computation is compromised

% Initialize next cells
cleanData = {[]};
badChan = {[]};

% Loop over the two subjects within the dyad
for subi = 1:2

    %refresh cfg
    cfg = [];
    cfg.method = 'summary'; %(?)
    %set limits for the graph
    cfg.alim = 5e-5;

    %For now, I want to store it! So I set the output cleanP1;
    %Base on literature the choice of he criterion (e.g. variance, kurtosis...
    %it's always specified in papers). Click and drag around point to mark for
    %rejection. After pressing 'Quit', selection is rejected.
    %after cleaning, you can type in command window to see the header of the
    %new structure
    cleanData{subi} = ft_rejectvisual(cfg,trialData{subi});

    % Write ch and trials down, plus indicator of interest e.g. variance; 
    % Acoording to FT tutorials, save in another file;
    % Anyway, save at least bad channels because you'll need to add back  
    % for comparing across subjects,e.g. with interpolation)
    % in my case, for having single trials, I'd rather remove bad channels; when
    % you have a lot of trials, you tend to remove trials
    badChan{subi} = trialData{subi}.label( find(~ismember(trialData{subi}.label, cleanData{subi}.label)) );
    %find a way to do the same for bad trials... alternative to write them down
    %by hand
    badTrial  = {[] ; []};

    %Rereference AFTER bad chan removal!!! So the eccessive noise does not
    %leak into the average; 
    cfg = [];
    cfg.reref='yes';
    cfg.refchannel={'all'}; 
    cleanData{subi} = ft_preprocessing(cfg , cleanData{subi});

    
end

clear contData;
clear filtData;

%I should sae trial rejected, also cause afterward i have to cancel that
%index from the trials of the partner



%% ICA
% This is good to submit to cluster looping over each subject. Then take
% one day to remove components from all subjects (visual inspection), then
% re-submit what comes next.

cfg = [];
% With 'runica', remove reference for being a linear combination
% of the sources; 
cfg.method       = 'runica';
cfg.channel      = {'all' '-CPz'};
cfg.trials       = 'all';
cfg.numcomponent = 'all';
cfg.demean       = 'no';

components = {[]};

tic;
parfor subi = 1:2 
    %note this calculation is not deterministic: it might sliiightly vary if
    %you iterate; anyway, it's not so important in practice
    components{subi} = ft_componentanalysis(cfg, cleanData{subi});

    %In principle you can continue analyzing the data on the component level by
    %using components as argument for ft_ functions instead of 'data'
   
end
clear trialData;

% Set path to save dataset
cd;
% Save dataset so far
save(['Coordination_' num2str(pair) '_ICAmatrix']);


% Check-point
close all;
clear all;



%% Artifact removal (by visual inspection)

close all;
clear all;


% ALWAYS check that you are working with the right set!!
pair = 2;

load(['Coordination_' num2str(pair) '_ICAmatrix']);

% In this case, analyses at electrode level: ICA for artifact removal
for subi = 1:2
  
    % plot the components topography; config inside the loop
    cfg = [];
    cfg.component = 1:round(length(components{subi}.label)/2);       % specify the component(s) that should be plotted
    cfg.layout    = layout; % specify the layout file that should be used for plotting (and just for plotting)
    cfg.comment   = 'no';
    ft_topoplotIC(cfg, components{subi})
    % further inspection of the time course of the components; oculomotor artifacts are particularly
    % evident 
    cfg = [];
    cfg.layout = layout; % specify the layout file that should be used for plotting
    cfg.viewmode = 'component';
    ft_databrowser(cfg, components{subi});

    % Wait for key press in command window before plotting subject2
    pause;

end

%% Remove the bad components and backproject the data

% Insert picked bad components
badComp = {[] ; []};
% Initialize cell for data pruned from artifactual components
postIca = {[]}; 

for subi = 1:2

    % Configuration within the loop, as different components across subjects 
    cfg = [];
    cfg.component = badComp{subi}; % to be removed component(s)
    postIca{subi} = ft_rejectcomponent(cfg, components{subi}, cleanData{subi});

end


%% Compare data before and after 

% Select channel and trial 
chan2plot = 'Fp1';         % frontal channels are ideal for assessing blink removal
chanIdx = find(strcmp(postIca{1}.label, chan2plot));
trialIdx = 18;             % pick a trial
exInterval = 1000:35000;   % pick a time window

figure;
for subi = 1:2
    subplot(2,1,subi);
    plot(cleanData{subi}.trial{trialIdx}(chanIdx,exInterval));
    hold on
    plot(postIca{subi}.trial{trialIdx}(chanIdx,exInterval) , 'r');
    legend('Before ICA removal' , 'After ICA removal');
    ylim([-150 150]);
    title(['Channel ' num2str(chan2plot)  ' - Subject ' num2str(subi) ' - Trial ' num2str(trialIdx)]);
end

%% Interpolation

% Preallocate for interpolated data
interpData = {[]};

%for exporting
eegData = {[]}; 
eegChan = {[]};


for subi = 1:2


    %prepare neighbours structure based on layout
    cfg = [];
    cfg.method = 'triangulation';
    cfg.layout = layout;
    cfg.feedback = 'yes';
    neighbours = ft_prepare_neighbours(cfg);    %load('elec1010_neighb.mat'); does not have M1 and M2


    %interpolation (it's important when it comes to compare across subjects)
    cfg = [];
    cfg.neighbours = neighbours; %defined in previous step
    cfg.method = 'weighted'; %default
    cfg.layout = layout; %use layout??
    cfg.badchannel = badChan{subi};
    %cfg.senstype = 'EEG';

    %NB: Mike x Cohen says interpolation should be done before ICA for optimizing
    %the algorithm ; in case you do it before, use cleanData as input
    %dataset (in general, mind wathever is the most updated); actually,
    %there are schools of thought... some people do this just fine and have
    %strong arguments on why this is better...

    % ... BUT interpolated channels are NOT independent anymore! so makes
    % more sense after ICA
    interpData{subi} = ft_channelrepair(cfg, postIca{subi});


    %EXPORT
    %pre-allocate
    eegData{subi} = zeros([size(interpData{1}.trial{1}) , ntrials]);
    %Assign, sorting conditions by the right order per each dyad
    for condi = 1:4
        condidx = find(condorder(pair,:)==condi);
    for triali = 1:10
        eegData{subi}(:,:,blocks(condidx,triali)) = interpData{subi}.trial{blocks(condi,triali)};
    end
    end
    eegTime = interpData{subi}.time(:,1);
    eegChan{subi} = interpData{subi}.label;
    srate = interpData{subi}.fsample;
    pnts = length(eegTime);


end
    
    % Export for integrated Matlab framework
    cd ; %set path
    save(['Coordination_' num2str(pair) '_eeg.mat'] , 'eegData' , 'eegTime' , 'eegChan' , 'srate');
    
    %store post-ICA data in fieldtrip format
    clear eegData;
    clear cleanData;  % or cleanData = {};
    clear postIca;
    save(['Coordination_' num2str(pair) '_ICApruned']);


    close all
    clear




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%Data are clean!





