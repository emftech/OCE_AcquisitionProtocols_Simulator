function T = AcqSim_SummaryTable(Metrics)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer (batch layer written with Claude Code)
% Date of creation (DD/MM/YY): 20-10-2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thesis-style summary table (rows = metrics, columns = estimators).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fn = fieldnames(Metrics);
MetricNames = [ ...
    "Mean Wave Phase Speed [m/s] ± std - Background (thesis disk)"
    "Accuracy Error Ea [%] - Background (thesis disk)"
    "Precision Error Ep [%] - Background (thesis disk)"
    "Mean Wave Phase Speed [m/s] ± std - Annulus (needle excluded)"
    "Accuracy Error Ea [%] - Annulus"
    "Precision Error Ep [%] - Annulus"
    "Coverage R2>=0.7 [%] - Annulus"
    "Mean Wave Phase Speed [m/s] ± std - Far field (r >= R + 1 mm, edge-free)"
    "Accuracy Error Ea [%] - Far field"
    "Precision Error Ep [%] - Far field"
    "Coverage R2>=0.7 [%] - Far field"];
T = table(MetricNames, 'VariableNames', {'Metric'});
for k = 1:numel(fn)
    M = Metrics.(fn{k});
    col = [ ...
        string(sprintf('%.4f ± %.4f', M.Bkg.mean, M.Bkg.std))
        string(sprintf('%.4f', M.AccuracyError.Ea_Bkg))
        string(sprintf('%.4f', M.PrecisionError.Ep_Bkg))
        string(sprintf('%.4f ± %.4f', M.Annulus.mean, M.Annulus.std))
        string(sprintf('%.4f', M.Annulus.Ea))
        string(sprintf('%.4f', M.Annulus.Ep))
        string(sprintf('%.1f', M.Annulus.Coverage))
        string(sprintf('%.4f ± %.4f', M.FarField.mean, M.FarField.std))
        string(sprintf('%.4f', M.FarField.Ea))
        string(sprintf('%.4f', M.FarField.Ep))
        string(sprintf('%.1f', M.FarField.Coverage))];
    T.(erase(fn{k}, ["MB_","SBM_"])) = col;
end
end
