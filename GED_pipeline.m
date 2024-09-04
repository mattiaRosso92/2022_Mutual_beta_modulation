% Generalized Eigendecomposition  

clear 
close all
clc

cd('/Volumes/G-DRIVE USB-C/eegData/Drifting_metronomes/chill_ICA');


% Design settings
ndyads  = 14; % dyads
nsubs   = 2;  % subjects per dyad
nconds  = 4;  % conditions
ncycles = 10;
nchans  = 64;
nsteps  = 64; % metronome's steps (i.e., cycles segments)
condlabels = {'Visual Coupling' , 'Visual Control' , 'Audio Coupling' , 'Audio Control'};

% Upload EEGlab structure for proper order of Chans (provides topoplotIndie() with input) 
load('chanlocs_waveguard64'); %
% Assign channel locations (as defined in standard_waveguard64.elc);
chanlocs_elc = {EEG.chanlocs.labels}';   

% ROIs selection
% chans2keep = chanlocs_elc;

chans2keep = {'T7','C3','Cz','C4','T8','CP5','CP1','CP2','CP6','P7','P3','Pz','P4',...
 'P8', 'POz','O1','Oz','O2','C5','C1','C2','C6','CP3','CPz','CP4','P5',...
 'P1','P2','P6','PO5','PO3','PO4','PO6','TP7','TP8','PO7','PO8'	}; % centro-posterior
                      
% chans2keep = {'Fp1','Fpz','Fp2','F7','F3' ,'Fz','F4','F8','FC5','FC1','FC2','FC6', ...
%     'AF7','AF3','AF4','F5','F1','F2','F6','FC3','FCz','FC4'}; % centro-frontal
%    
% Import behavioural data, for segmentation
load onsets_debounced
% ... and re-arrange in cell
tidx = cell(ndyads,2,nconds);
for dyadi = 1:ndyads  
    for condi = 1:nconds
        
            tidx{dyadi,1,condi} = onsets_debounced(dyadi).sub1{condi};
            tidx{dyadi,2,condi} = onsets_debounced(dyadi).sub2{condi};
            
        
    end
end 


%Initialize stuff (dyads,partners)
%components time series
compts    = cell(ndyads,nsubs,nconds); 
%and spectra
powr = cell(ndyads,nsubs,nconds); %raw
snr  = cell(ndyads,nsubs,nconds); %normalized signal-to-noise units

%covariance matrices (higher-level structure; cell content will be initialized later)
covS = cell(ndyads,nsubs,nconds);
covR = cell(ndyads,nsubs,nconds);
%and average
covS_avg = cell(ndyads,nsubs,nconds);
covR_avg = cell(ndyads,nsubs,nconds);

%GED
map   = zeros(ndyads,2,nconds,nchans); %forward model, for scalp topography; always full size! (64)
evals = zeros(ndyads,2,nconds,length(chans2keep)); %eigen values
evecs = zeros(ndyads,2,nconds,length(chans2keep),length(chans2keep)); %eigen vectors



%% GED parameters

% Center frequency
targetfrex = 20;
% Plateau filter
frange      = [targetfrex-2 targetfrex+2];   %frequency range
trans_width = .15; %percentage for transition zone (slope cut-off)
filt_ord    = 20; %filter order %NB I used 20 for beta band... why so?

% shrinkage proportion
shr = .01;

% time window for covarince matrix
twin = [-100 500]; % around tapping onsets


%% Apply RESS


for dyadi = 1:ndyads

    %Import and set
    load (['Coordination_' num2str(dyadi) '_eeg_chill']); %the order of conditions is already standardized 1234 for all subjects (maybe double-check in preprocessing script)
    time = eegTime{1};
    blockidx = [1:10 ; 11:20 ; 21:30 ; 31:40]; %smart indexing for dividing conditions

  
    % Synthetic empty channels for dyad 4
    if dyadi == 4
        for subi = 1:nsubs
        % Concatenate flat channels
        eegData{subi} =  cat(1,eegData{subi},zeros(length(chanlocs_elc)-length(eegChan{subi}),size(eegData{subi},2),size(eegData{subi},3)));       
        % Add missing channels
        eegChan{subi} = cat( 1,eegChan{subi},chanlocs_elc(~ismember(chanlocs_elc,eegChan{subi})) );
        end
    end
    
    % REORDER HERE!! (Already tested, but good to double-CHECK for analysis with Ole)
    for subi = 1:2
        [~ , chanidx] = ismember(chanlocs_elc , eegChan{subi});
        eegChan{subi} = eegChan{subi}(chanidx);     %re-sort channel labels
        eegData{subi} = eegData{subi}(chanidx,:,:); %re-sort data accordingly
    %A = [eegChan{1} , eegChan{1}(chanidx)]
    end
    
    % number of timepoints in filter
    pnts = length(time)*ncycles;
    
    % Frequency resolution (in Cohen's script, it was 0.05)
    frexres = 1/max(time); %Rayleigh frequency
    % FFT parameters
    nfft = ceil( srate/frexres ); %length sufficiently high to guarantee resolution at denominator
    %nfft = srate / (1/floor(max(time)*ntrials)); %Overkill??
    hz   = linspace(0,srate,nfft); %vector of frequencies; shortcut: all the frex above nyquist are not valid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Test RESS temporal filter parameters
    % filter parameters
%     frex = [ targetfrex-neig targetfrex targetfrex+neig ];
%     stds = [ fwhm_neig fwhm_targ fwhm_neig ]; %full width at half-maximum
%     % Visualize filters (function by Cohen)
%     RESSfilterFGxTest(pnts,srate,frex,stds,100);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    for subi = 1:2
                
            % Get indexes of channels to include in cov matrix
            post_logic = cellfun(@(c)strcmp(c,chanlocs_elc),chans2keep,'UniformOutput',false);
            post_idx = zeros(length(post_logic),1);
            for i = 1:length(post_logic)
                if find(post_logic{i}) >= 1 % correction ad-hoc for dyad4
                post_idx(i) = find(post_logic{i});
                end
            end
            post_idx = post_idx(post_idx~=0); % correction ad-hoc for dyad4
        
        for condi = 1:nconds % use blockidx to loop
            
            %if length(eegChan{subi}) == 64   %needs to match size of chanlocs.... SOLVE reprocessing when it's not the case: now I can interpolate anything        
            disp(['GED: dyad ' num2str(dyadi) '. Subject ' num2str(subi) '. Condition: ' condlabels{condi}])

            %temporary variable for time index
            tempt = tidx{dyadi,subi,condi};

            %initialize covariance matrix with 3rd dimension (based on N segments)
            covS{dyadi,subi,condi} = zeros(length(post_idx),length(post_idx),length(tempt));
            covR{dyadi,subi,condi} = zeros(length(post_idx),length(post_idx),length(tempt));
            %distance from grand-average (temporary variable)
            covS_dist = zeros(length(tempt),1);
            covR_dist = zeros(length(tempt),1);      

            %Assign trials selection to temporary input variable
            broad_long = eegData{subi}(post_idx,:,blockidx(condi,:)); %broadband signal
            broad_long = reshape(broad_long,length(post_idx),[],1);
            % Narrow band signal (choose your filter)
            narrow_long = filter_plateau (broad_long, srate, frange, trans_width, filt_ord,  1); %Plateau-shaped filter
            %narrow_long = filterFGx( broad_long , srate , targetfrex , fwhm_targ , 1); %Gaussian filter
                for ti = 1:length(tempt)

                    %Segment and assign time-window for covariance matrix
                    data = narrow_long(:,(twin(1)+tempt(ti)):(twin(end)+tempt(ti)-1));

%% S covariance matrix                    
                    % Compute covariance matrix
                    covS{dyadi,subi,condi}(:,:,ti) = cov(data');

%% R covariance matrices

                    % Assign broad-band segment to input data
                    data = broad_long(:,(twin(1)+tempt(ti)):(twin(end)+tempt(ti)-1));          
                    % ... and reshape                        
                    %data = reshape(data,length(eegChan{subi}),[],1);
                    % compute covariance
                    covR{dyadi,subi,condi}(:,:,ti) = cov(data');

                    % Apply shrinkage to covR (optional)
                    covR{dyadi,subi,condi}(:,:,ti) = (1-shr)*covR{dyadi,subi,condi}(:,:,ti)...
                        + shr*mean(eig(covR{dyadi,subi,condi}(:,:,ti)))*eye(size(covR{dyadi,subi,condi}(:,:,ti)));




                end


%% Grand-Average matrices

                %Compute average covariance matrices
                covS_avg{dyadi,subi,condi} = squeeze(mean(covS{dyadi,subi,condi} , 3)); %along 3D dimension
                covR_avg{dyadi,subi,condi} = squeeze(mean(covR{dyadi,subi,condi} , 3));

                %Remove 'bad' matrices based on distance
                for ti = 1:length(tempt)
                    %Compute distance
                    covS_dist(ti) = sqrt(trace(covS{dyadi,subi,condi}(:,:,ti)*covS_avg{dyadi,subi,condi}));
                    covR_dist(ti) = sqrt(trace(covR{dyadi,subi,condi}(:,:,ti)*covR_avg{dyadi,subi,condi}));

                end
                %Normalize distance (temporary variable)
                covS_z = (covS_dist - mean(covS_dist)) / std(covS_dist);
                covR_z = (covR_dist - mean(covR_dist)) / std(covR_dist);
                z_thresh = 2.3;  %~.01  

%                 %Visualize distances            
%                 figure(100+nconds+1) , clf
%                 subplot(211)
%                 plot(covS_z,'bs-','linew',2,'markerfacecolor','w')
%                 hold on
%                 plot(covR_z,'ks-','linew',2,'markerfacecolor','w')
%                 legend('CovS z-scores' , 'CovR z-scores')
%                 yline(z_thresh,'--');
%                 yline(-z_thresh,'--');
%                 title('Distances before removal')

                %find indexes outliers
                toofarS = abs(covS_z) > z_thresh;
                toofarR = abs(covR_z) > z_thresh;   

                %remove distances just for visualization, using logical indexing
                covS_z(toofarS) = NaN;
                covR_z(toofarR) = NaN;
%                 %visualize difference
%                 subplot(212)
%                 plot(covS_z,'bs-','linew',2,'markerfacecolor','w')
%                 hold on
%                 plot(covR_z,'ks-','linew',2,'markerfacecolor','w')
%                 legend('CovS z-scores' , 'CovR z-scores')
%                 yline(z_thresh,'--');
%                 yline(-z_thresh,'--');
%                 title('Distances after removal')

                %Remove actual outlier matrices, using logical indexing
                covS{dyadi,subi,condi}(:,:,toofarS) = NaN;
                covR{dyadi,subi,condi}(:,:,toofarR) = NaN;

                %Re-compute grand-averages (whatch out the NaNs)
                covS_avg{dyadi,subi,condi} = squeeze(mean(covS{dyadi,subi,condi} , 3 , 'omitnan')); %along 3D dimension
                covR_avg{dyadi,subi,condi} = squeeze(mean(covR{dyadi,subi,condi} , 3 , 'omitnan'));

                
%                 %Visualize cleaned Grand-Average covariance matrices
%                 figure(100+condi) , clf
%                 title(condlabels{condi})
%                 subplot(121)
%                 imagesc(covS_avg{dyadi,subi,condi})
%                 axis square
%                 title('Covariance Matrix S')
%                 subplot(122)
%                 imagesc(covR_avg{dyadi,subi,condi})
%                 axis square
%                 title('Covariance Matrix R')


%% Generalized Eigendecomposition (GED) 

                % Compute GED  
                [tempvecs,tempvals] = eig(covS_avg{dyadi,subi,condi},covR_avg{dyadi,subi,condi}); %assign to temporary vars 
                [tempvals,sidx]  = sort(diag(tempvals),'descend'); %sort components
                tempvecs = tempvecs(:,sidx);          %vectors of weights, for weighted average
                tempvals = 100*tempvals/sum(tempvals);   %express in percentage of explained variance (for eigenspectrum)

                % Pick the component by ranking of associated eigenvalue
                evalidx = 3;
               
                % Compute filter forward model and flip sign
                map(dyadi,subi,condi,post_idx) = tempvecs(:,evalidx)'*covS_avg{dyadi,subi,condi}; 
                [~,maxchan] = max(abs(map(dyadi,subi,condi,:)));
                map(dyadi,subi,condi,:) = map(dyadi,subi,condi,:)*sign(map(dyadi,subi,condi,maxchan));


                %% Component time series

                % Compute component time series (i.e., apply spatial filter to
                % "raw" data)
                compts{dyadi,subi,condi} = tempvecs(:,evalidx)'*reshape(broad_long,length(post_idx),[]); %WHY RESHAPE? CHECK THAT IT'S ALREADY THE RIGHT SHAPE
                compts{dyadi,subi,condi} = reshape(compts{dyadi,subi,condi},[],ncycles);

                % Power spectrum averaged over trials
                powr{dyadi,subi,condi} = mean( abs(fft(compts{dyadi,subi,condi},nfft,1)/pnts).^2 ,2);

                % SNR spectrum
%                 skipbins =  5; % .5 Hz
%                 numbins  = 20+skipbins; % 2 Hz
%                 hzdist = 2; %distance in hz on each side of centre frequency
%                 hzsize = 2; %size of reference windows in hx
%                 skipbins = round(frexres/hzdist); 
%                 numbins  = hzsize+skipbins; 
                %consider we are sliding the window over ALL the spectrum. Don't worry
                %about the alpha range; the important is that you have the
                %spike emerging from CLOSE neighbours over the whole
                %spectrum. Alpha emerges already due to GED
                skipbins = 20; %with full resolution, this is ~.5hz (same as Cohen) - interval hz(21)-hz(1); 
                numbins  = 20+skipbins; %%with full resolution, this is ~2hz
                snr{dyadi,subi,condi} = zeros(size(powr{dyadi,subi,condi}));


                for hzi=numbins+1:length(hz)-numbins-1

                    % SNR over all time points and conditions
                    numer = powr{dyadi,subi,condi}(hzi);
                    denom = mean(powr{dyadi,subi,condi}([hzi-numbins:hzi-skipbins hzi+skipbins:hzi+numbins]) );
                    snr{dyadi,subi,condi}(hzi) = numer./denom;

                end    

                
                %Adjust size of dyad4 concatenating 0s for missing electrodes
%                 
%                 if dyadi ==4 
%                    tempvals = [tempvals ; 0 ; 0 ; 0 ; 0 ; 0] ;
%                    tempvecs = [ [tempvecs ; repmat(zeros(1,length(tempvals)-5), [5 1])] , [repmat(zeros(length(tempvals),1), [1 5])] ] ;                  
%                 end
%                 
                
                % Store temporary variables
                evals(dyadi,subi,condi,:)   = tempvals;
                evecs(dyadi,subi,condi,:,:) = tempvecs;
                

            %end

                      
        clc

        end
        
    end

end

%clear EEG structure for topoplot and temp variable
clear data
clear EEG


%% Plot visual control 
% this visualization is meant for 'Stability Index' only; uncomment what
% follows when you want to analyze dyads in hyperscanning

% Upload EEGlab structure to provide topoplotIndie() with input 
load('chanlocs_waveguard64'); %

% Set condition: visual uncoupled (for 'Stability Index paper')
condition = 2;

for dyadi = 1:ndyads
    figure(dyadi), clf

    for subi = 1:2
        
        % Topography
        subplot(2,3,subi*3-2)
        topoplotIndie(map(dyadi,subi,condition,:),EEG.chanlocs,'numcontour',0,'electrodes','off');
        title('Component activation map') 
        % SNR spectrum
        subplot(2,3,subi*3-1)
        plot(hz,smooth(powr{dyadi,subi,condition},1),'k','linew',1.3)
        set(gca,'xlim',[10 40])     
        xlabel('Frequency (Hz)'), ylabel('Signal-to-noise ratio (%)')
        title('SNR spectrum')
        axis square
        sgtitle('Grand-average GED output (N = 28)')
        % Eigenspectrum
        subplot(2,3,subi*3);
        plot(squeeze(evals(dyadi, subi, condition,:)),'s-','markerfacecolor','k','linew',2)
        xlabel('Component'), ylabel('Eigenvalue')
        title('Eigenspectrum')
        axis square

    end

end


% %% Plotting 
% 
% %NB!! instead of difference, subplot 23X, with modality(2) and
% %mapCoupl,mapCtrl,SNR/pow
% 
% 
% smth = 1; %smoothing factor; don't oversmooth or you miss the whole point!!
% 
% % Initialize difference forward model, per modality
% vis_diff = zeros(nchans,ndyads,2);                       
% aud_diff = zeros(nchans,ndyads,2);  
% 
% for dyadi = 1:ndyads
%     figure(dyadi) , clf
% 
%     for subi = 1:2
% 
%         if sum(map(dyadi,subi,condi,:)) ~= 0   %check forward model validity (i.e., matching 64chans)
%             load('chanlocs_waveguard64');
%             %Trick: load an EEGLAB data structure to use
%             %locs as argument for topoIndie; clear EEG at the end
%             %NB when you interpolate chans, you append them at the end...
%             %someday, check indiePlot() and find out whether it is based on
%             %labels or order of channels. If it is the former, then I am
%             %completely fine now
% % 
% %             %To check quality of decomposition with eigenspectrum
% %             plot(squeeze(evals(pairi,subi,condi,1:20)),'s-','markerfacecolor','k','linew',2)
% %             xlabel('Component'), ylabel('% explained variance')
% %             title([ num2str(targetfrex) ' Hz' ])
% %             axis square
% 
% 
%             % Compute difference forward model, per modality
%             % NB: reflect on what is being mapped, and on the sign
%             vis_diff(:,dyadi,subi) = squeeze(map(dyadi,subi,1,:)) - squeeze(map(dyadi,subi,2,:));                       
%             aud_diff(:,dyadi,subi) = squeeze(map(dyadi,subi,3,:)) - squeeze(map(dyadi,subi,4,:));                       
%             
%             
%             % Plot forward models per modality (Coupling - Control)
%             subplot(2,4,subi*4-3)
%             topoplotIndie(squeeze(map(dyadi,subi,1,:)),EEG.chanlocs,'numcontour',0,'electrodes','off');
%             title([ num2str(targetfrex) ' Hz' ])
%             title('Visual coupling') 
%             
%             subplot(2,4,subi*4-2)
%             topoplotIndie(squeeze(map(dyadi,subi,3,:)),EEG.chanlocs,'numcontour',0,'electrodes','off');
%             title([ num2str(targetfrex) ' Hz' ])
%             title('Audio coupling') 
%            
%             
%             % Plot SNR per modality (Coupling - Control)
%             subplot(2,4,subi*4-1)
%             plot(hz,smooth(snr{dyadi,subi,1},smth),'k','linew',1.5) %vizcoupl
%             hold on
%             plot(hz,smooth(snr{dyadi,subi,2},smth),'r--','linew',1.2) %vizctrl
%             yline(1,'k--'); xline(targetfrex,'k--');
%             %set(gca,'xlim',[targetfrex-targetfrex targetfrex+targetfrex])     
%             xlim([0 5])
%             xlabel('Frequency (Hz)'), ylabel('SNR')
%             title(['Visual SNR: ' num2str(targetfrex) ' Hz' ])
%             %legend('Coupling', 'Control')
%             axis square
%             hold off
%             %To check raw spectrum
% %             plot(hz,powr{pairi,subi,1},'k','linew',1.5)
% %             hold on
% %             plot(hz,powr{pairi,subi,2},'r--','linew',1.2)
% %             set(gca,'xlim',[0 max(targetfrex+10)])
% %             xlabel('Frequency (Hz)'), ylabel('Raw power')
% %             title([ num2str(targetfrex) ' Hz' ])
%             %legend(['Coupling' 'Control'])
%             axis square
% 
% 
%             subplot(2,4,subi*4)
%             plot(hz,smooth(snr{dyadi,subi,3},smth),'k','linew',1.5) %audcoupl
%             hold on
%             plot(hz,smooth(snr{dyadi,subi,4},smth),'r--','linew',1.2) %audctrl
%             yline(1,'k--'); xline(targetfrex,'k--');
%             %set(gca,'xlim',[targetfrex-targetfrex targetfrex+targetfrex])     
%             xlim([0 5])
%             xlabel('Frequency (Hz)'), ylabel('SNR')
%             title(['Audio SNR: ' num2str(targetfrex) ' Hz' ])
%             %legend('Coupling', 'Control')
%             axis square
%             hold off
% %    %To check raw spectrum
% %             plot(hz,powr{pairi,subi,3},'k','linew',1.5)
% %             hold on
% %             plot(hz,powr{pairi,subi,4},'r--','linew',1.2)
% %             set(gca,'xlim',[0 max(targetfrex+10)])
% %             xlabel('Frequency (Hz)'), ylabel('Raw power')
% %             title([ num2str(targetfrex) ' Hz' ])
% %             %legend(['Coupling' 'Control'])
% %             axis square
% %             hold off
% %                    
%         end
%         
%     end
%         
% end
% 
% 
%% Grand averages (N=28)

ndyads = 14;
%ndyads = ndyads - 1;
targetfrex = 20;

%Initialize
map_gravg = zeros(nchans,nconds); 
powr_gravg = zeros(size(powr{1,1,1},1),nconds);
snr_gravg  = zeros(size(snr{1,1,1},1),nconds);
snr_single = zeros(length(hz),ndyads*2,nconds);
% Compute 
map_temp = squeeze(reshape(map,[],1,nconds,nchans)); %N=28, concatenate partners
for condi = 1:nconds
    map_gravg(:,condi) = squeeze(mean(map_temp(:,condi,:),1));
    powr_temp  = [powr{:,1,condi} powr{:,2,condi}]; %reshape cells into array with N=28 () concatenate partners
    powr_gravg(:,condi) = squeeze(mean(powr_temp,2)); %compute grand average
    snr_single(:,:,condi)  = [snr{:,1,condi} snr{:,2,condi}]; %reshape cells into array with N=28 () concatenate partners
    snr_gravg(:,condi) = squeeze(mean(snr_single(:,:,condi),2)); %compute grand average
end


figure(ndyads+1), clf
for condi = 1:nconds
    htopo = subplot(2,2,condi);
    topoplotIndie(map_gravg(:,condi),EEG.chanlocs,'numcontour',0,'electrodes','off');
    title(condlabels{condi}) 
    if targetfrex > 12 %if in Beta range
        caxis(htopo,[-.05 .05])
    else
        caxis(htopo,[-.1 .1])
    end
end

smth = 1; %smoothing factor

figure(ndyads+2), clf
%Raw power
subplot(221)
plot(hz,smooth(powr_gravg(:,1),smth),'k','linew',1.5)
hold on            
plot(hz,smooth(powr_gravg(:,2),smth),'r--','linew',1.2)
xline(targetfrex,'k--');
set(gca,'xlim',[0 50])
xlabel('Frequency (Hz)'), ylabel('Raw power')
title(['Visual raw power: ' num2str(targetfrex) ' Hz' ])
legend('Coupling' , 'Control')
hold off

subplot(222)
plot(hz,smooth(powr_gravg(:,3),smth),'k','linew',1.5)
hold on            
plot(hz,smooth(powr_gravg(:,4),smth),'r--','linew',1.2)
xline(targetfrex,'k--');
set(gca,'xlim',[0 50])
xlabel('Frequency (Hz)'), ylabel('Raw power')
title(['Audio raw power: ' num2str(targetfrex) ' Hz' ])
legend('Coupling' , 'Control')

%SNR
subplot(223)
plot(hz,smooth(snr_gravg(:,1),smth),'k','linew',1.5) %vizcoupl
hold on
plot(hz,smooth(snr_gravg(:,2),smth),'r--','linew',1.2) %vizctrl
yline(1,'k--'); xline(targetfrex,'k--');
set(gca,'xlim',[0 30])     
xlabel('Frequency (Hz)'), ylabel('SNR')
title(['Visual SNR: ' num2str(targetfrex) ' Hz' ])
legend('Coupling' , 'Control')
hold off

subplot(224)
plot(hz,smooth(snr_gravg(:,3),smth),'k','linew',1.5) %audcoupl
hold on
plot(hz,smooth(snr_gravg(:,4),smth),'r--','linew',1.2) %audctrl
yline(1,'k--'); xline(targetfrex,'k--');
set(gca,'xlim',[0 30])     
xlabel('Frequency (Hz)'), ylabel('SNR')
title(['Audio SNR: ' num2str(targetfrex) ' Hz' ])
legend('Coupling' , 'Control')
hold off

% Average eigenspectrum
evals28 = squeeze(reshape(evals(:,:,condition,:),ndyads*2,1,1,[]));
evalsAvg = mean(evals28,1);
