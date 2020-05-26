%% load in patient

load('lungPatient0.mat')

importOptions.inhExhDisc = false;
importOptions.numPhases = 5;

%% run this code, stop after bins are made

ct = matRad_binFrames2Phases(ct,importOptions);

%% now run this code to generate figure

figure
hold on
plot(t,3*(x_XCAT-x_XCAT(1)),'k')
plot([0 5],3*(repmat(lBounds,2,1)-x_XCAT(1)),'b')

xlabel('time / ss')
ylabel('tumour distance / mmmm')
grid on
set(gca,'YMinorTick','on','XMinorTick','on')
plot(t(ind_x_l_XCAT),3*(x_XCAT(ind_x_l_XCAT)-x_XCAT(1)),'kp','MarkerFaceColor','yellow')





figurePath  = 'C:\Users\eric\Carleton University\OneDrive - Carleton University\Carleton\PhD Project\Thesis\figures\generator\';
fname       = 'XCATTumourTrajectory_raw';
fullpath    = [figurePath fname '.tex'];

matlab2tikz('filename',fullpath,'interpretTickLabelsAsTex',true,'parseStrings',false,'noSize',true,'extraAxisOptions','minor y tick num=3,minor x tick num=3,width=\linewidth,height=\linewidth','showInfo',false);
