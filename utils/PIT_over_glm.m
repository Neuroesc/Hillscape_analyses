% function PIT_over_glm
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DESCRIPTION
% PIT_fig_1_v1  figure script written for:
% Grieves, Duvelle and Taube (202X) 
%
% See also: GIT_audit

% HISTORY:
% version 1.0.0, Release 06/11/23 Code conception
%
% Author: Roddy Grieves
% Dartmouth College, Moore Hall
% eMail: roddy.m.grieves@dartmouth.edu
% Copyright 2023 Roddy Grieves

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Heading 3
%% >>>>>>>>>>>>>>>>>>>> Heading 2
%% >>>>>>>>>> Heading 1
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> INPUT ARGUMENTS CHECK
    % Create figure
    fig_now = figure('Units','pixels','Position',[50 50 210.*3 297.*3],'visible','on');
    set(gcf,'InvertHardCopy','off'); % gives the figure a grey background but means it will save white lines as white    
    set(gcf,'color','w'); % makes the background colour white
    fs = [15 10];

x = 0:10:3000;
y = 0:10:1500;
F = 0.299;
wav = cos(2*pi*F*x+pi);
wav = ((wav+1)./2)*450;
[yy,xx] = ndgrid(y,x);
zz = repmat(wav,size(xx,1),1);
g = gradient(wav);
gg = repmat(g,size(xx,1),1);

clumaa.terrain_score = NaN(size(clumaa,1),1);
for mm = 1:size(clumaa,1)
    if clumaa.partn(mm)~=2
        continue
    end
    mnow = clumaa.ratemap_surficial{mm};
    gnow = imresize(gg,size(mnow),"nearest");

    r = corr(mnow(:),gnow(:),"Type","Pearson","Rows","pairwise");
    clumaa.terrain_score(mm,1) = r;
end

clumaa = PIT_correlations(config,pidx,clumaa,posdata,0);

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY

% Dependent variable
repetition_score = clumaa.repetition_score(pidx & clumaa.partn==2,1); % repetition score, all place cells, ridges

% Predictors
speed_score = clumaa.speed_score(pidx & clumaa.partn==2,1); % Nx1 vector, correlation values
azimuthal_directionality = clumaa.hd_3d_info(pidx & clumaa.partn==2,3); % Nx1 vector, spatial info content
left_vs_right_stability = clumaa.azimuth_stability(pidx & clumaa.partn==2,1);
pitch_directionality = clumaa.hd_3d_info(pidx & clumaa.partn==2,7); % Nx1 vector, spatial info content
up_vs_down_stability = clumaa.pitch_stability(pidx & clumaa.partn==2,1);
terrain_measure = abs(clumaa.terrain_score(pidx & clumaa.partn==2,1)); % Nx1 vector, terrain correlation


% figure
% scatter(repetition_score,speed_score,30,'k','filled','o');
% refline
% keyboard

    % Pack predictors and dependent variable into a table
    T = table(speed_score(:), azimuthal_directionality(:), left_vs_right_stability(:), pitch_directionality(:), up_vs_down_stability(:),terrain_measure(:), repetition_score(:),'VariableNames', {'speed_score','azimuthal_directionality','left_vs_right_stability','pitch_directionality','up_vs_down_stability','terrain_measure','repetition_score'});

    % Full GLM (linear regression)
    mdl = fitglm(T,'repetition_score ~ speed_score + azimuthal_directionality + left_vs_right_stability + pitch_directionality + up_vs_down_stability + terrain_measure','Distribution', 'normal', 'Link', 'identity');
    plotPartialRegression(fig_now,mdl);

    % Display outputs
    disp('--- Full Model Summary ---');
    disp(mdl);


% Get full model R²
R2_full = mdl.Rsquared.Ordinary;

% Partial R² loop (from previous code)
predictors = mdl.PredictorNames;
partialR2 = zeros(numel(predictors),1);

for i = 1:numel(predictors)
    reducedPredictors = setdiff(predictors, predictors{i});
    formula = sprintf('repetition_score ~ %s', strjoin(reducedPredictors, ' + '));
    mdl_reduced = fitglm(mdl.Variables, formula);
    
    SSE_reduced = mdl_reduced.SSE;
    SSE_full = mdl.SSE;
    
    % Partial R² (unique contribution of predictor i)
    partialR2(i) = (SSE_reduced - SSE_full) / SSE_reduced;
end

% % of total R² explained
percentR2 = (partialR2 / R2_full) * 100;

% Results table
T_contrib = table(predictors, partialR2, percentR2,'VariableNames', {'Predictor','PartialR2','PercentOfTotalR2'})

% keyboard
%% >>>>>>>>>> Save the overall figure
    if 1
        fname = [config.fig_dir '\Fig S11a.png']; 
        if fast_figs
            frame = getframe(gcf); % fig is the figure handle to save
            [raster, raster_map] = frame2im(frame); % raster is the rasterized image, raster_map is the colormap
            if isempty(raster_map)
                imwrite(raster, fname);
            else
                imwrite(raster, raster_map, fname); % fig_file is the path to the image
            end            
        else
            exportgraphics(gcf,fname,'BackgroundColor',[1 1 1],'Colorspace','rgb','Resolution',res);  
        end
        close(gcf);   
    end


 









function plotPartialRegression(fig,mdl)
    % mdl: fitglm or fitlm object
    
    X = mdl.Variables;           % full table of predictors + dependent
    y = mdl.Variables{:, mdl.ResponseName}; % dependent variable
    predictors = mdl.PredictorNames;        % list of predictor names
    predictor_plot_names = {'Speed score','Directionality','Directional stability','Tilt stability','Tilt directionality','Terrain score'};
    coefTable = mdl.Coefficients;              % GLM coefficient table

    figure(fig);
    tiledlayout(4,3)
    nPred = numel(predictors);
    
    for i = 1:nPred
        % Current predictor
        x = X.(predictors{i});
        
        % Regress y on all other predictors (excluding current one)
        otherPreds = setdiff(predictors, predictors{i});
        tbl_y = X(:, otherPreds);
        lm_y = fitlm(tbl_y, y);
        x_resid = lm_y.Residuals.Raw;
        
        % Regress current predictor on all others
        lm_x = fitlm(tbl_y, x);
        y_resid = lm_x.Residuals.Raw;
        
        % Plot residuals
        nexttile
        % subplot(ceil(sqrt(nPred)), ceil(sqrt(nPred)), i);
        alpha = 0.5;
        scatter(x_resid, y_resid, 20,'k','filled','MarkerFaceAlpha',alpha,'MarkerEdgeColor','none');
        ylabel([predictor_plot_names{i} ' (residuals)'],'Interpreter',"none");
        if i>3
            xlabel(['Repetition score (residuals)'],'Interpreter',"none");
        else
            xlabel('  ')
        end
        % title(['Partial regression: ' predictors{i}]);
        hold on;
        
        % Fit line to residuals
        ax = gca;
        ax.XLim = [-0.5 1];
        ax.YLim(2) = ax.YLim(2)*1.15;
        p = polyfit(x_resid, y_resid, 1);
        x_fit = linspace(ax.XLim(1), ax.XLim(2), 100);
        y_fit = polyval(p, x_fit);
        plot(x_fit, y_fit, 'r-', 'LineWidth', 2);

        % Extract GLM results for this predictor
        row = coefTable(strcmp(coefTable.Row, predictors{i}), :);
        beta = row.Estimate;
        se = row.SE;
        tval = row.tStat;
        pval = row.pValue;
        
        if pval < .05
            % Scientific notation: mantissa × 10^{exponent}
            [mantissa, exponent] = sprintf('%.2e', pval);  % format in scientific
            parts = regexp(sprintf('%.2e', pval), '([-+]?\d*\.\d+)e([-+]?\d+)', 'tokens');
            mantissa = parts{1}{1};
            exponent = str2double(parts{1}{2});
            pStr = sprintf('\\it{p} = %s \\times 10^{%d}', mantissa, exponent);
        else
            % Normal decimal, strip leading zero
            pStr = sprintf('\\it{p} = %.3f', pval);
            pStr = regexprep(pStr, '0\.', '.');
        end
        
        % Annotation string
        txt = sprintf('\\beta = %.3f, SE = %.3f\nt = %.2f, %s', beta, se, tval, pStr);
        
        % Place annotation at the very top of the plot
        text(0.5, 1, txt, 'HorizontalAlignment', 'center', 'VerticalAlignment','middle','FontSize', 7, 'Interpreter','tex','Units','normalized');

    end
end



