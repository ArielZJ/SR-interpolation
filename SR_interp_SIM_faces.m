%% SIMULATE SIMILARITY RATING DATA AND USE INTERPOLATION METHOD

%Generate simulated data that contains TWO categories (here, male or female
%faces) of which NO EXACT stimulus pairings are ever shown.

% matSizes = The size of the original stimulus matrix (used to generate
%simulated data) in the form matSizes*matSizes
% possDiff = Set the possible SR (within category)
% levelsOfVar = Noise amongst all SR (between category)

%% Collect participant data, find group average similarity matrices and plot
clear all; close all;

matSizes = 4; % 8; 16; %set size of original SR matrix
maxDissimVal = 7; %maximum dissimilarity value
possDiff = [0,7; 0.5,6.5; 1,6; 1.5,5.5; 2,5; 2.5,4.5; 3,4; 3.25,3.75; 3.5, 3.5]; %max SR values for category 1 and 2
levelsOfVar = [2.5, 2, 1.5, 1, 0.5, 0.1, 0.05]; %level of noise amongst ALL SR

%possDiff = [0,7; 2,5; 3.5, 3.5];
%levelsOfVar = [2.5, 1,0.05];

countPlot = 1; %idx
triangle = 0; %set whether to calculate triangle inequality or not (1 = yes)

for mm = 1:length(matSizes)
    
    matSize = matSizes(mm);
    
    for cc = 1:length(possDiff)
        
        meanDiff(cc) = abs(possDiff(cc,1)-possDiff(cc,2));
        simVal = possDiff(cc,1);
        diffVal = possDiff(cc,2);
        
        for kk = 1:length(levelsOfVar)
            
            useVar = levelsOfVar(kk)
            
            savename = ['SIMULATED_' num2str(matSize) 'x' num2str(matSize) '_Matrix_w_MeanSRDiff_' ...
                num2str(meanDiff(cc)) '_Var_' num2str(useVar) '.mat']
            
            %Generate variance/noise for simulated SR data
            n = matSize*matSize;
            varRange = [0 useVar];
            genData = rand(n,1)*range(varRange)+min(varRange);
            base_SR_sim = reshape(genData, [matSize,matSize]);
            
            %Generate the max similar and dissimilar values (based on mean difference)
            withinSR_SIM = repmat(simVal,[matSize/2,matSize/2]); %set max SR sim value
            withinSR_DIFF= repmat(diffVal,[matSize/2,matSize/2]); %set max SR diff value
            allSRs = [withinSR_SIM, withinSR_DIFF;
                withinSR_DIFF', withinSR_SIM']; %combine same/diff SRs
            
            %Combine mean SRs and noise/variance
            simulated_Data_AttPeriph = allSRs+base_SR_sim;
            
           %Show the results from the simulation
            figure(1111);
            subplot(size(possDiff,1),size(levelsOfVar,2),countPlot)
            imagesc(simulated_Data_AttPeriph,[0 7.5]); colorbar;
            
            %% Setup empty matrices to be filled
            %Create empty matrix for 16 x 16 comparisons
            finalDissim_AttPeriph = NaN(matSize*2,matSize*2);
            
            %Create empty new matrices to input to final matrix at the end
            dissim_AttPeriph_TopLeftQuandrant = NaN(matSize/2,matSize/2);
            dissim_AttPeriph_TopRightQuandrant = NaN(matSize/2,matSize/2);
            dissim_AttPeriph_BottomLeftQuandrant = NaN(matSize/2,matSize/2);
            dissim_AttPeriph_BottomRightQuandrant = NaN(matSize/2,matSize/2);
            
            dissim_AttPeriph_TopRight =  simulated_Data_AttPeriph;
            
            % Gather matrices - FIRST we take the transpose to give the bottom left corner
            dissim_AttPeriph_BottomLeft = transpose(dissim_AttPeriph_TopRight);
            
            %Add the above to the final dissimilarity matrix
            finalDissim_AttPeriph(1:matSize,matSize+1:matSize+matSize) = dissim_AttPeriph_TopRight;
            finalDissim_AttPeriph(matSize+1:matSize+matSize,1:matSize) = dissim_AttPeriph_BottomLeft;
            
            %% FIRST COMPILE SRs FOR THE OVERALL BOTTOM RIGHT QUADRANT (OF 16 X 16 MATRIX)
            % That is, compile an 8 x 8 matrix that will fit into bottom right quadrant
            % Gather data for TOP LEFT this 8 x 8 quadrant
            
            colComb = nchoosek(1:matSize/2,2); %get pairs to compare based on matrix size (i.e. half)
            
            %FOR THE TOP LEFT OF THE MATRIX
            for ii = 1:size(colComb,1)
                
                vRA = dissim_AttPeriph_TopRight(1:matSize/2,colComb(ii,1)); %select vector RA
                vRB = dissim_AttPeriph_TopRight(1:matSize/2,colComb(ii,2)); %select vecotr RB
                
                %Check stdev
                svRA = std(vRA); svRB = std(vRB);
                if svRA == 0; vRA(2) = vRA(2) + 0.00001; end %if corr = 0 = err. Add small amount of noise.
                if svRB == 0; vRB(2) = vRB(2) + 0.00001; end
                
                corr_vLA_vLB = corr(vRA, vRB); %compute correlation
                SRating = -1 * (corr_vLA_vLB-1) / 2 * maxDissimVal; %scale SR to match corr
                
                addSR = colComb(ii,:); %coords to input data
                addSRb = flip(addSR); %same coords (symmetry) for A/B = B/A
                
                dissim_AttPeriph_TopLeftQuandrant(addSR(1),addSR(2)) = SRating; %add rating
                dissim_AttPeriph_TopLeftQuandrant(addSRb(1),addSRb(2)) = SRating; %add rating
            end
            
            %Gather data for  BOTTOM RIGHT of the 8 x 8 matrix
            colComb = nchoosek(1:matSize/2,2)+matSize/2; % assign columns to compare (as above)
            intoMatPos = nchoosek(1:matSize/2,2); %add this value, to this position in new quadrant
            
            for ii = 1:size(colComb,1)
                vRA = dissim_AttPeriph_TopRight(matSize/2+1:matSize/2+matSize/2,colComb(ii,1)); %extract vector RA
                vRB = dissim_AttPeriph_TopRight(matSize/2+1:matSize/2+matSize/2,colComb(ii,2)); %extract vector RB
                
                %Check stdev
                svRA = std(vRA); svRB = std(vRB);
                if svRA == 0; vRA(2) = vRA(2) + 0.00001; end
                if svRB == 0; vRB(2) = vRB(2) + 0.00001; end
                
                corr_vLA_vLB = corr(vRA, vRB);
                SRating = -1 * (corr_vLA_vLB-1) / 2 * maxDissimVal;
                
                addSR = intoMatPos(ii,:);
                addSRb = flip(addSR); %coords to input data
                
                dissim_AttPeriph_BottomRightQuandrant(addSR(1),addSR(2)) = SRating;
                dissim_AttPeriph_BottomRightQuandrant(addSRb(1),addSRb(2)) = SRating;
                
            end
            
            % Gather data for BOTTOM LEFT of 8 x 8 matrix
            possCoords = [1:matSize/2]; %select these coordinates
            
            [A1,A2] = ndgrid(possCoords);
            colComb = [A2(:),A1(:)];
            
            %Assign the colums to compare
            colComb(:,2) = colComb(:,2)+matSize/2;
            
            %Assign the rows/cols in the new matrix to place these values
            intoMatPos = [A2(:),A1(:)];
            
            for ii = 1:size(colComb,1)
                
                vRA = dissim_AttPeriph_TopRight(1:matSize/2,colComb(ii,1));
                vRB = dissim_AttPeriph_TopRight(1:matSize/2,colComb(ii,2));
                
                %Check stdev
                svRA = std(vRA); svRB = std(vRB);
                if svRA == 0; vRA(2) = vRA(2) + 0.00001; end
                if svRB == 0; vRB(2) = vRB(2) + 0.00001; end
                
                corr_vLA_vLB = corr(vRA, vRB);
                SRating = -1 * (corr_vLA_vLB-1) / 2 * maxDissimVal;
                
                addSR = intoMatPos(ii,:);
                
                dissim_AttPeriph_BottomLeftQuandrant(addSR(1),addSR(2)) = SRating;
                
            end
            
            %TOP RIGHT quandrant will be transpose of BOTTOM LEFT
            dissim_AttPeriph_TopRightQuandrant = transpose(dissim_AttPeriph_BottomLeftQuandrant);
            
            
            %% NOW GATHER EACH QUADRANT FOR A FULL 8 x 8 matrix
            % THIS REPRESENTS THE FINAL BOTTOM RIGHT OF THE NEW MATRIX!
            final_Dissim_newBOTTOMRIGHT = [dissim_AttPeriph_TopLeftQuandrant, dissim_AttPeriph_TopRightQuandrant; ...
                dissim_AttPeriph_BottomLeftQuandrant, dissim_AttPeriph_BottomRightQuandrant];
            
            finalDissim_AttPeriph(matSize+1:matSize+matSize,matSize+1:matSize+matSize) = final_Dissim_newBOTTOMRIGHT;
            
            
            %% NOW COMPILE SRs FOR THE OVERALL TOP LEFT QUADRANT (OF THE FINAL 16 X 16 MATRIX)
            % That is, compile an 8 x 8 matrix that will fit into TOP LEFT quadrant
            
            %Create empty new matrices to input to final matrix at the end
            dissim_AttPeriph_TOP_TopLeftQuandrant = NaN(matSize/2,matSize/2);
            dissim_AttPeriph_TOP_TopRightQuandrant = NaN(matSize/2,matSize/2);
            dissim_AttPeriph_TOP_BottomLeftQuandrant = NaN(matSize/2,matSize/2);
            dissim_AttPeriph_TOP_BottomRightQuandrant = NaN(matSize/2,matSize/2);
            
            % Gather data for TOP LEFT this 8 x 8 quadrant
            rowComb = nchoosek(1:matSize/2,2);
            
            %FOR THE TOP LEFT OF THE MATRIX
            for ii = 1:size(rowComb,1)
                
                vRA = dissim_AttPeriph_TopRight(rowComb(ii,1),1:matSize/2)';
                vRB = dissim_AttPeriph_TopRight(rowComb(ii,2),1:matSize/2)';
                
                %Check stdev
                svRA = std(vRA); svRB = std(vRB);
                if svRA == 0; vRA(2) = vRA(2) + 0.00001; end
                if svRB == 0; vRB(2) = vRB(2) + 0.00001; end
                
                corr_vLA_vLB = corr(vRA, vRB);
                SRating = -1 * (corr_vLA_vLB-1) / 2 * maxDissimVal;
                
                addSR = rowComb(ii,:);
                addSRb = flip(addSR); %coords to input data
                
                dissim_AttPeriph_TOP_TopLeftQuandrant(addSR(1),addSR(2)) = SRating;
                dissim_AttPeriph_TOP_TopLeftQuandrant(addSRb(1),addSRb(2)) = SRating;
            end
            
            %Gather data for  BOTTOM RIGHT of the 8 x 8 matrix
            rowComb = nchoosek(1:matSize/2,2)+matSize/2;
            intoMatPos = nchoosek(1:matSize/2,2);
            
            for ii = 1:size(rowComb,1)
                
                vRA = dissim_AttPeriph_TopRight(rowComb(ii,1),1:matSize/2)';
                vRB = dissim_AttPeriph_TopRight(rowComb(ii,2),1:matSize/2)';
                
                %Check stdev
                svRA = std(vRA); svRB = std(vRB);
                if svRA == 0; vRA(2) = vRA(2) + 0.00001; end
                if svRB == 0; vRB(2) = vRB(2) + 0.00001; end
                
                corr_vLA_vLB = corr(vRA, vRB);
                SRating = -1 * (corr_vLA_vLB-1) / 2 * maxDissimVal;
                
                addSR = intoMatPos(ii,:);
                addSRb = flip(addSR); %coords to input data
                
                dissim_AttPeriph_TOP_BottomRightQuandrant(addSR(1),addSR(2)) = SRating;
                dissim_AttPeriph_TOP_BottomRightQuandrant(addSRb(1),addSRb(2)) = SRating;
            end
            
            % Gather data for BOTTOM LEFT of 8 x 8 matrix
            possCoords = [1:matSize/2];
            
            [A1,A2] = ndgrid(possCoords);
            rowComb = [A2(:),A1(:)];
            
            %Assign the rows to compare
            rowComb(:,2) = rowComb(:,2)+matSize/2;
            
            %Assign the rows/cols in the new matrix to place these values
            intoMatPos = [A2(:),A1(:)];
            
            for ii = 1:size(rowComb,1)
                
                vRA = dissim_AttPeriph_TopRight(rowComb(ii,1),1:matSize/2)';
                vRB = dissim_AttPeriph_TopRight(rowComb(ii,2),1:matSize/2)';
                
                %Check stdev
                svRA = std(vRA); svRB = std(vRB);
                if svRA == 0; vRA(2) = vRA(2) + 0.00001; end
                if svRB == 0; vRB(2) = vRB(2) + 0.00001; end
                
                corr_vLA_vLB = corr(vRA, vRB)
                SRating = -1 * (corr_vLA_vLB-1) / 2 * maxDissimVal;
                
                addSR = intoMatPos(ii,:)
                
                dissim_AttPeriph_TOP_BottomLeftQuandrant(addSR(1),addSR(2)) = SRating;
                
            end
            
            %TOP RIGHT quandrant will be transpose of BOTTOM LEFT
            dissim_AttPeriph_TOP_TopRightQuandrant = transpose(dissim_AttPeriph_TOP_BottomLeftQuandrant);
            
            %% NOW GATHER EACH QUADRANT FOR A FULL 8 x 8 matrix
            % FINAL QUADRANT TO ADD INTO THE COMPLETE 16 x 16 matrix
            final_Dissim_new_TOPLEFT = [dissim_AttPeriph_TOP_TopLeftQuandrant, dissim_AttPeriph_TOP_TopRightQuandrant; ...
                dissim_AttPeriph_TOP_BottomLeftQuandrant, dissim_AttPeriph_TOP_BottomRightQuandrant];
            
            finalDissim_AttPeriph(1:matSize,1:matSize) = final_Dissim_new_TOPLEFT;
            
            figure(2)
            subplot(size(possDiff,1),size(levelsOfVar,2),countPlot)
            imagesc(finalDissim_AttPeriph,[0 7]); colorbar;
            
           countPlot = countPlot+1;


        end
        
    end
end

