%% load results

clear all
load stats

nSubPerPosPhase_opt = 1;
nSubPerVelPhase_opt = 2;
nTimeFracs_opt      = 6;
r_opt               = 9;

rArr = round(26./fResampleArr);

path = 'C:\Users\eric\Carleton University\OneDrive - Carleton University\Carleton\PhD Project\Thesis\figures\generator\';

%% nTimeFracs vs r !!!!!!

figure

nTimeFracsVec  = 1:10;
rVec           = [26 9 4 3 2 1];

[nTimeFracsVec_grid,rVec_grid]  = meshgrid(nTimeFracsVec,rVec);

pSumImage = zeros(size(nTimeFracsVec_grid));

for i = 1:numel(rVec)
    for j = 1:numel(nTimeFracsVec)
        pSumImage(i,j) = pSumArr(nSubPerPosPhaseArr == nSubPerPosPhase_opt & nSubPerVelPhaseArr == nSubPerVelPhase_opt & ...
            nTimeFracsArr == nTimeFracsVec_grid(i,j) & rArr == rVec_grid(i,j));
    end
end

imagesc(pSumImage)
xticklabels(cellstr(num2str(nTimeFracsVec')))
yticklabels(cellstr(num2str(rVec')))
xlabel('N_T')
ylabel('r')
colorbar
caxis([0 0.7])

fname = 'pSum_nTimeFracs_vs_r_raw';
fullpath = [path fname '.tex'];
matlab2tikz('filename',fullpath,'interpretTickLabelsAsTex',true,'parseStrings',false,'noSize',true,'extraAxisOptions','minor y tick num=4,minor x tick num=4,width=\linewidth,height=\linewidth','showInfo',false);

%% nTimeFracs vs nSubPerPosPhase

figure

nTimeFracsVec       = 1:10;
nSubPerPosPhaseVec  = fliplr(1:10);

[nTimeFracsVec_grid,nSubPerPosPhaseVec_grid]  = meshgrid(nTimeFracsVec,nSubPerPosPhaseVec);

pSumImage = zeros(size(nTimeFracsVec_grid));

for i = 1:numel(nSubPerPosPhaseVec)
    for j = 1:numel(nTimeFracsVec)
        pSumImage(i,j) = pSumArr(nSubPerPosPhaseArr == nSubPerPosPhaseVec_grid(i,j) & nSubPerVelPhaseArr == nSubPerVelPhase_opt & ...
            nTimeFracsArr == nTimeFracsVec_grid(i,j) & rArr == r_opt);
    end
end

imagesc(pSumImage)
xticklabels(cellstr(num2str(nTimeFracsVec')))
yticklabels(cellstr(num2str(nSubPerPosPhaseVec')))
xlabel('N_T')
ylabel('N_x')
colorbar
caxis([0 0.7])

fname = 'pSum_nTimeFracs_vs_nSubPerPosPhase_raw';
fullpath = [path fname '.tex'];
matlab2tikz('filename',fullpath,'interpretTickLabelsAsTex',true,'parseStrings',false,'noSize',true,'extraAxisOptions','minor y tick num=4,minor x tick num=4,width=\linewidth,height=\linewidth','showInfo',false);


%% nTimeFracs vs nSubPerVelPhase !!!

figure

nTimeFracsVec       = 1:10;
nSubPerVelPhaseVec  = fliplr(1:10);

[nTimeFracsVec_grid,nSubPerVelPhaseVec_grid]  = meshgrid(nTimeFracsVec,nSubPerVelPhaseVec);

pSumImage = zeros(size(nTimeFracsVec_grid));

for i = 1:numel(nSubPerVelPhaseVec)
    for j = 1:numel(nTimeFracsVec)
        pSumImage(i,j) = pSumArr(nSubPerPosPhaseArr == nSubPerPosPhase_opt & nSubPerVelPhaseArr == nSubPerVelPhaseVec_grid(i,j) & ...
            nTimeFracsArr == nTimeFracsVec_grid(i,j) & rArr == r_opt);
    end
end

imagesc(pSumImage)
xticklabels(cellstr(num2str(nTimeFracsVec')))
yticklabels(cellstr(num2str(nSubPerVelPhaseVec')))
xlabel('N_T')
ylabel('N_v')
colorbar
caxis([0 0.7])

fname = 'pSum_nTimeFracs_vs_nSubPerVelPhase_raw';
fullpath = [path fname '.tex'];
matlab2tikz('filename',fullpath,'interpretTickLabelsAsTex',true,'parseStrings',false,'noSize',true,'extraAxisOptions','minor y tick num=4,minor x tick num=4,width=\linewidth,height=\linewidth','showInfo',false);


%% nSubPerPosPhase vs r !!!

figure

nSubPerPosPhase = 1:10;
rVec            = [26 9 4 3 2 1];

[nSubPerPosPhase_grid,rVec_grid]  = meshgrid(nSubPerPosPhase,rVec);

pSumImage = zeros(size(nSubPerPosPhase_grid));

for i = 1:numel(rVec)
    for j = 1:numel(nSubPerPosPhase)
        pSumImage(i,j) = pSumArr(nSubPerPosPhaseArr == nSubPerPosPhase_grid(i,j) & nSubPerVelPhaseArr == nSubPerVelPhase_opt & ...
            nTimeFracsArr == nTimeFracs_opt & rArr == rVec_grid(i,j));
    end
end

imagesc(pSumImage)
xticklabels(cellstr(num2str(nSubPerPosPhase')))
yticklabels(cellstr(num2str(rVec')))
xlabel('N_x')
ylabel('r')
colorbar
caxis([0 0.7])

fname = 'pSum_nSubPerPosPhase_vs_r_raw';
fullpath = [path fname '.tex'];
matlab2tikz('filename',fullpath,'interpretTickLabelsAsTex',true,'parseStrings',false,'noSize',true,'extraAxisOptions','minor y tick num=4,minor x tick num=4,width=\linewidth,height=\linewidth','showInfo',false);


%% nSubPerVelPhase vs r

figure

nSubPerVelPhaseVec  = 1:10;
rVec                = [26 9 4 3 2 1];

[nSubPerVelPhaseVec_grid,rVec_grid]  = meshgrid(nSubPerVelPhaseVec,rVec);

pSumImage = zeros(size(nSubPerVelPhaseVec_grid));

for i = 1:numel(rVec)
    for j = 1:numel(nSubPerVelPhaseVec)
        pSumImage(i,j) = pSumArr(nSubPerPosPhaseArr == nSubPerPosPhase_opt & nSubPerVelPhaseArr == nSubPerVelPhaseVec_grid(i,j) & ...
            nTimeFracsArr == nTimeFracs_opt & rArr == rVec_grid(i,j));
    end
end

imagesc(pSumImage)
xticklabels(cellstr(num2str(nSubPerVelPhaseVec')))
yticklabels(cellstr(num2str(rVec')))
xlabel('N_v')
ylabel('r')
colorbar
caxis([0 0.7])

fname = 'pSum_nSubPerVelPhase_vs_r_raw';
fullpath = [path fname '.tex'];
matlab2tikz('filename',fullpath,'interpretTickLabelsAsTex',true,'parseStrings',false,'noSize',true,'extraAxisOptions','minor y tick num=4,minor x tick num=4,width=\linewidth,height=\linewidth','showInfo',false);


%% nSubPerVelPhase vs nSubPerPosPhase

figure

nSubPerVelPhaseVec  = 1:10;
nSubPerPosPhaseVec  = fliplr(1:10);

[nSubPerVelPhaseVec_grid,nSubPerPosPhaseVec_grid]  = meshgrid(nSubPerVelPhaseVec,nSubPerPosPhaseVec);

pSumImage = zeros(size(nSubPerVelPhaseVec_grid));

for i = 1:numel(nSubPerPosPhaseVec)
    for j = 1:numel(nSubPerVelPhaseVec)
        pSumImage(i,j) = pSumArr(nSubPerPosPhaseArr == nSubPerPosPhaseVec_grid(i,j) & nSubPerVelPhaseArr == nSubPerVelPhaseVec_grid(i,j) & ...
            nTimeFracsArr == nTimeFracs_opt & rArr == r_opt);
    end
end

imagesc(pSumImage)
xticklabels(cellstr(num2str(nSubPerVelPhaseVec')))
yticklabels(cellstr(num2str(nSubPerPosPhaseVec')))
xlabel('N_v')
ylabel('N_x')
colorbar
caxis([0 0.7])

fname = 'pSum_nSubPerVelPhase_vs_nSubPerPosPhase_raw';
fullpath = [path fname '.tex'];
matlab2tikz('filename',fullpath,'interpretTickLabelsAsTex',true,'parseStrings',false,'noSize',true,'extraAxisOptions','minor y tick num=4,minor x tick num=4,width=\linewidth,height=\linewidth','showInfo',false);
