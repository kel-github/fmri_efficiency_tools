%%% tools to compare efficiency of event related designs in fMRI 
%%% written by K. Garner, 2019 (c). Free to use and share, please cite.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% to run the below functions, first simulate the temporal onsets of your trials
% and use these onsets to specify a design matrix in SPM (e.g. SPM -> Batch
% -> SPM -> stats -> fmri model specification (design only)

clear all

%% add paths for functionality
addpath('~/Documents/spm12/'); % change to your spm path

%% define variables for the below functions
TRs = 1.5; % a vector of possible TR durations
iters = 1; % number of iterations
figName = 'TR%d_sim%d'; % root of all figure names
figLoc  = pwd; % folder in which to save figs
outPath = 'sim_spms/TR1500_TRsc1_sim1/'; % where the predefined spm file is saved

%% plot the regressors in time and freq domain with hi-pass cutoff and save
% use plot regressors function 
for iTRs = 1:numel(TRs)
    
   for iIters = 1:iters
       
      fName = sprintf(figName, TRs(iTRs)*1000, iIters);
      spmLoc  = outPath;
      load([spmLoc '/' 'SPM.mat'])
      spm = SPM;
      vocal = 1;
      plot_regressors( spm, figLoc, fName, vocal );      
   end
end

%% define contrast vectors
% for each regressor, I take the total number of regressors, the matrix
% denoting the design for that subset of regressors (e.g. a 4 x 3 design), 
% and the index of the subset of regressors used in the specified design.
% I then zero pad using the index, and save to a structure
tReg = 17; % total number of experimental regressors (excluding session regressors)
designs = { [4, 2], 2, 3, 2 }; % each element is the levels of each factor to be included in that design
c_idx = {[1, 2, 4, 5, 7, 8, 10, 11],...
          [13, 14],...
          [15, 16, 17], ...
          [16, 17]}; % each element denotes the regressors of of interest for each contrast design
conts = make_contrasts_KG( tReg, designs, c_idx );


%% compute efficiency
all_effs = cell(numel(TRs), iters);
for iTRs = 1:numel(TRs)
   TR = TRs( iTRs ); 
   
   for iIters = 1:iters   
     
        % load SPM file
       spmLoc  = outPath;
       load( [ spmLoc '/' 'SPM.mat' ] );
       all_effs{ iTRs, iIters } = run_efficiency_comp_KG(SPM, conts);
             
   end  
end

