%mainFolder contains all the cells of interest. each cell should have three
%channels, corresponding to the gefDomain, KRAS domain, and membrane
%marker. 
%savePath is the path where plots are saved
%mySegParameters if the path of the segmentationFile.

% mainFolder = '/endosome/archive/bioinformatics/Danuser_lab/shared/microdomain/sourceData/tiam1_DecemberVersion/'
% savePath = '/endosome/archive/bioinformatics/Danuser_lab/shared/microdomain/results/Year26Results_March/03_26/labelsGabrielAnalysisPDF/'
% myLabels = readmatrix('/endosome/archive/bioinformatics/Danuser_lab/shared/microdomain/sourceData/tiamLabels_2.csv')

mainFolder = ''
savePath = ''
mySegParameters = ''
cd(mainFolder)
myDir = dir(mainFolder);

mkdir(savePath)

selectedFolder = mySegParameters(:,1);
selectedMovies = mySegParameters(:,2);
selectedCell = mySegParameters(:,3);
selectedFrame = mySegParameters(:,4);
selectedThresh = mySegParameters(:,5);
myRatio = repmat(nan, [1,length(selectedMovies)])

%% Compute the ratio of activity inside/outide protrusions, and plot. 
for ff = 1 :length(selectedMovies)
    if (selectedCell(ff) < 100)
        folderIndex = selectedFolder(ff) + 2;
        wFolder = append(myDir(folderIndex).folder, '/', myDir(folderIndex).name);

        mySubfolder = dir(append(wFolder, '/field', string(selectedMovies(ff))))

        fileIndex = find(endsWith({mySubfolder.name}, 'during.tif'));
        toRead = append(mySubfolder(fileIndex).folder, '/', mySubfolder(fileIndex).name);

        gefDomain = readtiff(append(toRead, '_ch0.tif'));
        krasDomain = readtiff(append(toRead, '_ch1.tif'));
        membraneMarker = readtiff(append(toRead,'_ch2.tif'));

        for t = [2, selectedFrame(ff)] %only segment frames used in computation.
            myFrame = bwlabeln(bwareafilt(imfill(gefDomain(:,:,t) > selectedThresh(ff), 'holes'), [400, inf]), 4);

            myVector = myFrame(myFrame>0);
            tallies = tabulate(myVector);
            talliesSort = sortrows(tallies, 2, 'descend');
            % segCell(:,:,t) = (myFrame == mode(myFrame(myFrame>0)));
            segCell(:,:,t) = (myFrame == talliesSort(selectedCell(ff), 1));

        end
      
        krasNormalized = (10000 *krasDomain) ./ membraneMarker;

        krasNormalized_Protrusion = uint16(krasNormalized(:,:,selectedFrame(ff))) .* (uint16((segCell(:,:,2) == 0))) .* uint16((segCell(:,:,selectedFrame(ff)) == 1));
        krasNormalized_Outside = uint16(krasNormalized(:,:,selectedFrame(ff))) .* (uint16((segCell(:,:,2) == 1))) .* uint16((segCell(:,:,selectedFrame(ff)) == 1));

        krasNormalized_Protrusion(krasNormalized_Protrusion == 0) = nan;
        krasNormalized_Outside(krasNormalized_Outside ==0 ) = nan;

        valInside(ff) = mean(krasNormalized_Protrusion(krasNormalized_Protrusion>0), 'omitmissing'); %inside prootrusion
        valOutside(ff) = mean(krasNormalized_Outside(krasNormalized_Outside>0), 'omitmissing'); % outside protrusion

        % final ratio computation
        myRatio(ff) = valInside(ff)/valOutside(ff);
        
        % Plotting Section 
        winterMod = winter(100)
        winterMod(1,:) = [0,0,0]   
        maxValInside = ceil(max(krasNormalized_Protrusion(:))/50) * 50
        maxValOutside = ceil(max(krasNormalized_Outside(:))/50) * 50
        maxToPlot = max(maxValInside, maxValOutside);

        f = figure()
        tiledlayout(1,2)
        nexttile()
        imshow(krasNormalized_Protrusion, [0, maxToPlot])
        colormap(winterMod)
        colorbar()
        title(append('Protrusive Parts','Ratio',string(valInside(ff)/valOutside(ff))))

        nexttile()
        imshow(krasNormalized_Outside, [0, maxToPlot])
        colormap(winterMod)
        colorbar
        title('Non Protrusive Parts')

        saveas(f,append(savePath, '/sampleOutput',string(ff),'.pdf'))

    end

    close all
end


myTableOut = [valOutside(valOutside>0 .* (valInside>0)); valInside(valOutside >0 .* (valInside>0))]' 
writetable(uint16(myTableOut), 'rr_forGG_RasTable_02_08_26.txt')


