function [] = plot_regressors(spm, outDir, outNm, vocal)
% written by K. Garner, 2019, free to share, please cite and use
% responsibly
% adapted from spm function spm_DesRep()
% spmDir = SPM output structure
% outDir = location of where to save output figure
% outFM  = output filename stem
% vocal  = do you want the standard spm (design -> explore) plots for each session and 
%          regressor? (1). Or just the prob Spectra given the regressor plots
%          (0). 
SPM = spm;
%plot_regressors
if any(vocal)
    figure; % start new figure
end

Pspectra = zeros( 5, size( SPM.xBF.bf, 2), length(SPM.Sess(1).Fc ), length( SPM.Sess ) ); % for a bounds (hz) x parameters x regressors x session matrix

for s = 1:length(SPM.Sess)
    
    for i = 1:length(SPM.Sess(s).Fc)
        
        
        %-Trial-specific regressors - time domain
        %--------------------------------------------------------------------------
        Sess  = SPM.Sess;
        sX    = SPM.xX.X(Sess(s).row,Sess(s).col); % SPM.xX.X = design matrix, Sess(s).row = scan indices, Sess(s).col = effect indices for sessions
        rX    = sX(:,Sess(s).Fc(i).i); % rX = the design matrix for the session and regressors for the HRF under under question
        if any(vocal)
            subplot(2,2,1)
            plot(Sess(s).row,rX)
            xlabel('scan')
            ylabel('regressor[s]')
            title({'Time domain',['Regressors for ' Sess(s).Fc(i).name]})
            grid on
            axis tight
        end
        
        %-Trial-specific regressors - frequency domain
        %--------------------------------------------------------------------------
        gX    = abs(fft(rX)).^2; % fourier transform of rX for the session and regressors for the HRF under under question
        gX    = gX*diag(1./sum(gX)); % take 1./sum(gx) of each column and put on a diagonal matrix, multiply each value by that proportion - i.e. get p(x)
        q     = size(gX,1);
        Hz    = [0:(q - 1)]/(q*SPM.xY.RT); % hz in sample
        q     = 2:fix(q/2); % index of Hz up to the sampling frequency
        HPF   = SPM.xX.K(s).HParam; % get high pass
        
        % The below uses measures from the canonical HRF which has a typical FWHM of 6 seconds and duration of 32
        % seconds i.e. +/- 1/(32-3) & 1/(32+3)
        HFWHM = 6*.5;
        idx = dsearchn(Hz', [1/HPF 1/(32+HFWHM), 1/(32-HFWHM) .1]');
        Pspectra(:, :, i, s) =  [sum(gX(1:idx(1), :)); ...
                                 sum(gX(idx(1)+1:idx(2), :)); ...   
                                 sum(gX(idx(2)+1:idx(3), :));...
                                 sum(gX(idx(3)+1:idx(4),:));...
                                 sum(gX(idx(4)+1:end,:))];
     
        if any(vocal)
            subplot(2,2,2) 
            plot(Hz(q),gX(q,:)) % plot hz against frequency spectra of regressor
                    patch([0 1 1 0]/HPF,[0 0 1 1]*max(max(gX)),[1 1 1]*.9,'facealpha',.5);
            xlabel('Frequency (Hz)')
            ylabel('relative spectral density')
            h=title(['Frequency domain',sprintf('\n'), ' {\bf',num2str(HPF),'}', ...
                ' second High-pass filter'],'Interpreter','Tex');
            grid on
            axis tight
        end
        
        if any(vocal)
            % if trial (as opposed to trial x trial interaction)
            %--------------------------------------------------------------------------
            if length(Sess(s).U) >= i
                
                % Basis set and peristimulus sampling
                %----------------------------------------------------------------------
                subplot(2,2,3)
                dt   = Sess(s).U(i).dt;
                RT   = SPM.xY.RT;
                t    = [1:size(SPM.xBF.bf,1)]*dt;
                pst  = Sess(s).U(i).pst;
                plot(t,SPM.xBF.bf,pst,0*pst,'.','MarkerSize',16)
                str  = sprintf('TR = %0.2fs',RT);
                xlabel({'time {secs}' str sprintf('%0.0fms time bins',1000*dt)})
                title({'Basis set and peristimulus sampling' SPM.xBF.name})
                axis tight
                grid on
                
                % if a paramteric variate is specified
                %----------------------------------------------------------------------
                for p = 1:length(Sess(s).U(i).P)
                    
                    if Sess(s).U(i).P(p).h
                        
                        % onsets and parametric modulation
                        %------------------------------------------------------------------
                        subplot(2,2,4)
                        ons = Sess(s).U(i).ons;
                        plot(ons,Sess(s).U(i).P(p).P,'.','MarkerSize',8)
                        xlabel('time {secs}')
                        title('Parameters')
                        grid on
                        hold on
                        
                    end
                end
            end
        end
        
        %%
        if any(vocal)
            h = gcf;
            print(h, sprintf([outDir, '/' outNm '_reg_%s.pdf'], Sess(s).Fc(i).name), '-dpdf');
            clf(h, 'reset')
        end
    end
end

% now do stacked bar chart of frequency spectra proportions for exah sess x
% regressor
figure;
count_plot = 0;
for s = 1:length(Sess)
    for i = 1:size(Pspectra,2)
        count_plot = count_plot + 1;
        c = categorical({Sess(s).Fc.name});
        subplot(length(Sess), size(Pspectra,2), count_plot )
        h = bar(c,squeeze(Pspectra(:,i,:,s))', 'stacked');
        axis tight
    end
end
legend(h, {'<HP', 'HP-HFW', 'PEAK', 'HP-p1', 'REM' }, 'Location', 'Best');
h = gcf;
print(h, sprintf([outDir, '/' outNm 'pSpectra.pdf']), '-dpdf');
clf(h, 'reset')


