function [ conts ] = make_contrasts_KG( tReg, designs, cIdx )
% written by K. Garner, 2019, free to share, cite and use responsibly
% tReg = total number of regressors
% designs = cell, each element is a design vector where each i indicates
% how many levels in factor i, should be ordered according to spm
% assumptions,
% cIdx = cell, each element contains the row idxs which should be included
% in a given contrast

% conts = a cell where each element is a structure of the relevant
% contrasts to the corresponding elements in design
conts = {};

% first, generate contrasts of interest
% a contrast, e.g. a 4 x 2 design, or a 1 factor with 2 levels design etc
for iD = 1:length(designs)
    
    c(iD).c = spm_make_contrasts(designs{iD});
end

% for each c, make a matrix of zeros with reg columns and c rows, fill the
% appropriate columns of c with the c matrix
for iD = 1:length(c)
    fixCs = []; % empty variable in which the fixed contrasts for each design will be collated
    tmp = c(iD);
    for iC = 1:length(tmp.c)
        tmpC     = tmp.c(iC);
        contrast = tmpC.c;
        cMat     = zeros( size( contrast, 1 ), tReg );
        cMat( :, cIdx{ iD } ) = contrast;
        name     = tmp.c(iC).name;
        fixCs(iC).c    = cMat;
        fixCs(iC).name = name;
    end
    conts{ iD } = fixCs;
    
end