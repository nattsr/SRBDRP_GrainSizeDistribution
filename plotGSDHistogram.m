function [] = plotGSDHistogram(binCenter, grainSizePDF)
%plotGSDHistogram plot Grain Size Distribution output from computeGSDHistogram
%   Input Arguments
%   - binCenter     : a (49*1) vector, bin center in mm as specified in the
%                     code
%   - grainSizePDF  : a (49*1) vector, probability distribution function
%                     (PDF) - volume fraction

%   Revision 1: Feb 2018 Nattavadee Srisutthiyakorn



%% Program
% Plot
maxY = 0.5;
redColor = [0.6350 0.0780 0.1840];
bar(log10(binCenter), grainSizePDF,0.8, 'k');
% Set plot specification
xTickVec = -3:1;
%set(gca, 'XDir', 'reverse','Xtick', xTickVec, 'Xticklabel', 10.^xTickVec);
set(gca,'Xtick', xTickVec, 'Xticklabel', 10.^xTickVec);
ylabel(gca, 'Volume Fraction');
xlabel('Particle Diameter (mm)');
ylim([0 maxY])

box on;
hold on

% Plot Cutoff
cutoffNum = [2 1 0.5 0.25 0.125 0.0625 0.0039 0];
cutoffTxt = {'Gravel','VC Sand','C Sand','M Sand','F Sand','VF Sand','Silt','Clay'};
cutoffTxtLoc = [2.4 1.2 0.6 0.3 0.15 0.075 0.0048 0.003];
nTxt = length(cutoffNum);

for iTxt = 1:nTxt
    plot([log10(cutoffNum(iTxt)) log10(cutoffNum(iTxt))],[0 1],'--', 'color', [0.6350 0.0780 0.1840])
    h = text(log10(cutoffTxtLoc(iTxt)),0.25, cutoffTxt{iTxt});
    set(h, 'rotation', 90); h(1).Color = [0.6350 0.0780 0.1840];
end


end

