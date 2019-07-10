function [effs] = run_efficiency_comp_KG(spm, conts)
% K. Garner, 2019, free to share and use, please cite and share responsibly
% --------------------------------------------------------------------------
% this function takes the information from the spm structure (defined using
% SPM functionality - see batch_specify_design_only) and the pre-defined contrasts 
% (conts), and computes the design efficiency according to 
% http://imaging.mrc-cbu.cam.ac.uk/imaging/DesignEfficiency

% inputs: spm   = an SPM structure 
%         conts = a cell structure {}, each entry in the cell is a structure 
%                 with fields:
%         conts{1}.c   - Contrast matrix
%         conts{1}.name - Name
%
% output:
% effs = a cell array the same size as conts, with an efficiency score in each
% element, corresponding to the analagous contrast in conts

SPM   = spm; 
X     = SPM.xX.X; % get the design matrix
sess  = length(SPM.Sess); % n sessions to be used to concatenate the contrast vector for the linear algebra functions below
nHRFB = size(SPM.xBF.bf, 2); % n of basis functions for HRF

effs = cell(1, size( conts, 2 ) ); % number of efficiencies to compute

for iDes = 1:length( conts ) % scroll through contrasts for different
    
    % experimental design
    c_vec = conts{ iDes };  
    eff   = zeros(1, size( c_vec, 2) );
    for i = 1:length( c_vec )
        
        tmpC = c_vec(i).c;
        effs{iDes}(i) = compute_efficiency_KG( X, tmpC, sess, nHRFB );
    end
    
end