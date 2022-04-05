%% Function to do the anchoring analysis based on vessels
% 2021/11/12-
% Yunhao Bai & Bokai Zhu

% Features: 
   %'_cellRatio': instead of simple cell number encountered, use the cell overlap pixel ratio to reflect the cell distribution
   % use additional filter to solve the overlapping issue when adjacent tumor patches expand and overlap
   % use totalIon (closed) as mask to avoid boundaries that goes out of imaging area
   % solve the erosion issue which leads to elimination of some regions

%% main function
% -40 will etch through almost all vessels for ExPRESSO samples
distMin = -40; 
distMax = 20;

% import the panel of all markers
massDS = dataset('File','info/0421panel_ExMIBI-FFPE_brain_forBetty.csv','Delimiter',',');
% load the cell phenotypes
cellAnnotation = dataset('File','info/20211124_Brain_Expand_all_scaleSize_anno2.csv','Delimiter',',');
    
% path to the core folder
corePath = {'Sample_data/MFG_ExPRESSO_3x3_400_1024x1024_1ms_2depth'};
% folder of the subfolder for non-nuclear segmentation results
ezSegPath = {'ezSegResults_Betty_MFG_3x3'};
% folder of the cell segment results
segmentPath = {'H3CD450.3_mpp0.15_Mesmer'};  
coreNum = length(corePath);

% initiate the dataframe, if more than 1 core
%allDataL = [];

% process each core
for p=1:coreNum
    disp(['Currently: ',corePath{p}]);

    % initiate the dataframes for current core
    coreDataL = [];
    channelNum = length(massDS.Label);
    
    % load the ezSegment objects (other than Vessel): 
    % we do not need to distinguish between different fragments, so just import a simple mask
    t = imread([corePath{p},'/',ezSegPath{p},'/masks/Point1/Astrocytes2.tif']); %mind the .tif and .tiff
    Astrocytes2Map = double(t); %a 0 or 255 image
    t = imread([corePath{p},'/',ezSegPath{p},'/masks/Point1/Microglia2.tif']); %mind the .tif and .tiff
    Microglia2Map = double(t); %a 0 or 255 image 
    
    
    % load tiffs to recreate countsNoNoise
    imX = size(Astrocytes2Map,1); imY = size(Astrocytes2Map,2);
    % initiate the countsNoNoise, channelNum + 2x exSegment + 1x cellMask
    countsNoNoise = zeros(imX,imY,channelNum+2+1); 
    for i=1:channelNum
        t = imread([corePath{p},'/',massDS.Label{i}, '.tiff']); %mind the .tif and .tiff
        d = double(t);
        % imshow(d)
        countsNoNoise(:,:,i) = d;
    end
    % append the Astrocytes2Map and Microglia2Map after that
    countsNoNoise(:,:,end-2) = Astrocytes2Map;
    countsNoNoise(:,:,end-1) = Microglia2Map;

    % create an image+ region to avoid out-of-boarder regions, use totalIon
    % as the mask
    totalIon = countsNoNoise(:,:,1:end-3);
    totalIonIm = sum(totalIon,3);
    totalIonIm(totalIonIm == 0) = 0;
    SE = strel('square',256);
    % imclose to fill local holes
    totalIonIm2 = imclose(totalIonIm,SE);
    
    % load the cell mask 
    load([corePath{p},'/',segmentPath{p},'/segmentationParams.mat']);
    %get newLmod and labelIdentityNew
    labelIdentityNew2 = labelIdentityNew([1:end-1]); % fix bug resulting from previous script

    currentCellAnno = cellAnnotation(cellAnnotation.point_id == p,:); 
    %so the newLmod we get from above lines will have their phenotype anno
    countsNoNoise(:,:,end) = newLmod;
    
    % reshape for easier signal extraction and processing
    countsReshape = reshape(countsNoNoise,imX*imY,channelNum+2+1);
    
    % load the Vessel2 masks, etc
    load([corePath{p},'/',ezSegPath{p},'/objects_points/Point1/Vessel2_objData.mat']);
    
    % close the vessel if necessary, sometimes the vessel feature segment
    % has holes inside
    imCloseSet = {0}; %imclose images with these numbers
    SE = strel('diamond',imCloseSet{p});
    mapped_obj_ids2 = imclose(mapped_obj_ids,SE);
    objectNum = max(max(mapped_obj_ids2));
    % get the vesselSize, construct the vector
    vesselSTAT = regionprops(mapped_obj_ids2,'Area');
    vesselSizes = zeros(objectNum,1);
    for i=1:objectNum
        vesselSizes(i) = vesselSTAT(i).Area;
    end
    
    % treat each region individually to avoid overlapping during expansion
    %make individual vessel masks
    tempRegion = cell(1,objectNum);
    for i=1:objectNum
        tempRegionMatrix = mapped_obj_ids2;
        tempRegionMatrix(tempRegionMatrix ~= i) = 0;
        tempRegion{i} = tempRegionMatrix;
    end
    
    for dist = distMin:distMax
        disp(dist);
        % define the strcuture elements
        if dist ~= 0
            SE1 = strel('diamond',abs(dist));
            SE2 = strel('diamond',abs(dist)-1);
        end
        %stepRegionImAll = zeros(imX,imY);
        
        if dist == 0 %skip this situation
            continue;
        elseif dist < 0 %erosion
            % mind that with small distMin, the region can be eliminated
            % and generate NaN in follwing steps
            updatedIm = imerode(mapped_obj_ids2,SE1);
            updatedIm2 = imerode(mapped_obj_ids2,SE2);
            stepIm = updatedIm2 - updatedIm;
            updatedSTAT = regionprops(stepIm,'Area','PixelIdxList');
            % to avoid suddenly kill of the updatedSTAT when the length is
            % shorter than objectNum due to elimination of some tumors, match 
            % the length to avoid issues
            if length(updatedSTAT) < objectNum
                for v = length(updatedSTAT)+1:objectNum
                    updatedSTAT(v).Area = 0;
                    updatedSTAT(v).PixelIdxList = [];
                end
            end
        else %expansion
            % if expand too much, vessel may crossover with each others
            updatedSTAT = struct('Area',cell(1,objectNum),'PixelIdxList',cell(1,objectNum));
            for v = 1:objectNum
                updatedRegionIm = imdilate(tempRegion{v},SE1);
                updatedRegionIm2 = imdilate(tempRegion{v},SE2);
                stepRegionIm = updatedRegionIm - updatedRegionIm2;
                % filter to remove possible overlap with orginal other
                % vessels (not likely), only during expansion
                stepRegionIm(mapped_obj_ids2 ~= 0) = 0;
                % filter to remove vessel that expand out of the meaningful
                % region, only during expansion
                stepRegionIm(totalIonIm2 == 0) = 0;
                
                %stepRegionImAll = stepRegionImAll + v.*stepRegionIm; 
                % cannot simply add together, as will generated weird dots 
                % on overlapping pixels of adjacent regions
                tempSTAT = regionprops(stepRegionIm,'Area','PixelIdxList');
                updatedSTAT(v).Area = tempSTAT(v).Area;
                updatedSTAT(v).PixelIdxList = tempSTAT(v).PixelIdxList;
            end
        end
        
        % finally, extract the difference
        % vessel ID(1) + step ID(1) + vesselSize(1) + features/markers(labelNum) + ezSegment(2) + cells(4)
        data = zeros(objectNum,channelNum);
        dataScaleSize = zeros(objectNum,channelNum);
        stepSizes = zeros(objectNum,1);
        
        % coreID
        coreID = p*ones(objectNum,1);
        % vessel ID(1) + step ID(1)
        vesselID = [1:objectNum]; %remember to transpose
        stepID = dist*ones(objectNum,1); %remember to transpose
        
        % features/markers(labelNum)
        % for each label extract information
        for i=1:objectNum
            currData = countsReshape(updatedSTAT(i).PixelIdxList,1:end-3);
            data(i,1:channelNum) = sum(currData,1);
            dataScaleSize(i,1:channelNum) = sum(currData,1)/updatedSTAT(i).Area;
            stepSizes(i) = updatedSTAT(i).Area;
        end
        
        % Phenotypes(4)
        cellVector = zeros(objectNum,4);
        for i=1:objectNum
            currData = countsReshape(updatedSTAT(i).PixelIdxList,end);
            % get the cell ID
            cellList = unique(currData);
            cellList = cellList(cellList>0);
            if ~isempty(cellList)
                for j = 1:length(cellList)
                    [row, col] = find(currData == cellList(j));
                    currentCellNum = find(currentCellAnno.obj_id == cellList(j));
                    cellVector(i,currentCellAnno.nuMeta(currentCellNum)) = cellVector(i,currentCellAnno.nuMeta(currentCellNum)) + length(row)/stepSizes(i);
                end
            end
        end
        
%%
        % the fragment ratios, ezSegment(2)
        ezVector = zeros(objectNum,2);
        for i=1:objectNum
            currData = countsReshape(updatedSTAT(i).PixelIdxList,end-2:end-1);
            % Astrocytes2Map, then Microglia2Map 
            ezVector(i,:) = sum(currData,1)./(255*updatedSTAT(i).Area); %mind the maps are 0 or 255
        end
        
        % assemble the function 
        stepDataL = [coreID,vesselID',vesselSizes,stepID,stepSizes,dataScaleSize,cellVector,ezVector];
        % remove the row with Area==0 or have NaN value because of too
        % small distMin which leads to elimination of some regions
        stepDataL(any(isnan(stepDataL),2), :) = []; %here choose NaN
        % finally assemble the dataL for this core
        coreDataL = [coreDataL;stepDataL];
    end
    % if more than 1 core, assemble all the dataframes
    %allDataL = [allDataL;coreDataL];
    csvwrite(['vessel2_stepwise_cellRatio_P',num2str(p),'_',num2str(distMin),'to',num2str(distMax),'.csv'],coreDataL);
end

% if have multiple cores, write everything into a matrix
%csvwrite(['20211205_vessel2_stepwise_cellRatio',num2str(distMin),'to',num2str(distMax),'_boundaryExpand_.csv'],allDataL);




 




