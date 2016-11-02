%Thomas An
%A11257870
%ECE174 Coding Assigment 1 LPC Lossy Compression of Speech Signals
%Due 11/01/16
%Professor Ken Kreutz-Delgado

y = audioread('An_Thomas_A11257870_Project1.wav');
%sound(y);
%%
%STEP 2: Quantizing the Y signal directly and then storing the MSE, I went
%back and added the for loops to store multiple MSE data with varying R and
%Alpha values
%initializing variables to store blocks and samples
blocksize = 160;
totalcells = floor(length(y)/blocksize);
%getting and eliminating truncated data in y
remainderdata = mod(length(y),blocksize);
yremainder = y(1:(length(y)-remainderdata));
ydyn_range = max(yremainder) - min(yremainder);
y_length = length(yremainder);
%creating cell array to store everything and truncating remainders

%C is the cell array to store the signal
C = cell(totalcells,1);
%These two are to be used during quantization to retain original values
Cthresh = cell(totalcells,1);
Cthresh2 = cell(totalcells,1);
%k and j used to iterate through y
k = 1;                                
j = 160;
Ymean = cell(totalcells,1);
Ystddev = cell(totalcells,1);

%Breaking down signal into the cellarray
for itr=1:totalcells    
    %For loop storing audio file into 160 sample size cells
    C{itr,1} = y(k:j); 
    Cthresh{itr,1} = y(k:j);
    Cthresh2{itr,1} = y(k:j);
    %increment to move to next 160 elements in array
    k = k+160;
    j = j+160;
    
end

yq = cell(totalcells,1);
%MSE = (yremainder-yqcolumn)'*(yremainder-yqcolumn)/y_length;
MSEMatrixYQ = zeros(8,7);
%In this loop, we are quantizing the cell array with the Y stored in,
%and then storing the MSE into an MSEMatrix which will be used later
%for graphing purposes.
%We iterate through r = 1:8 and alpha = 0.5:3.5 because we want to gather
%multiple MSE values with carry r and alpha to be able to pick the right
%data points
for r = 1:1:8
    alphaindex = 1;
    L = 2^r;
    for alpha = 0.5:0.5:3.5
        yq = cell(totalcells,1);
        Cthresh = Cthresh2;
        for itr=1:totalcells
            Ymean{itr,1} = mean(Cthresh2{itr,1});
            Ystddev{itr,1} = std(Cthresh2{itr,1});
            %setting upper and lower bounds 1 stddev away
            %For loop eliminating outliers
            upperbound2 = Ymean{itr,1} + alpha*Ystddev{itr,1};
            lowerbound2 = Ymean{itr,1} - alpha*Ystddev{itr,1};
            for t=1:blocksize
                if Cthresh{itr,1}(t,1) > upperbound2
                    Cthresh{itr,1}(t,1) = upperbound2;
                elseif Cthresh{itr,1}(t,1) < lowerbound2
                    Cthresh{itr,1}(t,1) = lowerbound2;
                end
            end
        end
        %quantization of the Y
        yq = cell(totalcells,1);
        for i=1:totalcells
            q = (max(Cthresh{i,1})-min(Cthresh{i,1}))/(L-1);
            yq = round(Cthresh{i,1}/q)*q;
            Cthresh{i,1} = yq;
        end
        yqcolumn = cell2mat(Cthresh);
        %sound(yqcolumn);
        MSEloop = (yremainder-yqcolumn).'*(yremainder-yqcolumn)/y_length;
        MSEMatrixYQ(r,alphaindex) = MSEloop;
        alphaindex = alphaindex + 1;
    end
end

%%
%STEP 3: Creating and finding the toeplitz matrix and filter coefficient
%matrix
l1=10;
%cells for storing fillcoeffs and toepstores
toepstore = cell(totalcells, 1);
fillstore = cell(totalcells,1);
%setting initial condition for previous 

previous = zeros(blocksize, 1);
%this for loop iterate through previous block as well as the current block
%to create the column and rows for the toeplitz command
%the toeplitz matrices are then stored as well as the filter coefficients
for itr = 1:totalcells
    current = Cthresh{itr,1};
    col = [previous(end);current(1:end-1)];
    row = previous(end:-1:end-(l1-1));
    toepmatrix1 = toeplitz(col,row);
    fillstore{itr,1} = toepmatrix1\current;
    %fillstore{itr,1} = fillcoeff;
    toepstore{itr,1} = toepmatrix1;
    previous = current;
    if itr ~= totalcells
        current = C{itr+1,1};
    end
end
%%
%STEP 4: Finding residuals using the toepliz matrix and filter coefficients
%create e to store residuals of every block
e1 = cell(totalcells,1);
for itr = 1:totalcells
    e1{itr,1} = C{itr,1} - toepstore{itr,1}*fillstore{itr,1};
end
e1column = cell2mat(e1);
e1dyn_range = max(e1column) - min(e1column);
%%
%STEP 5: Find the residuals directly from the equation listed, then
%comparing the MSE values between that found in step 4
%initializing variables and matrices to store data and residuals
y_hat = cell(totalcells,1);
e2 = cell(totalcells,1);
%the next for loop does the summation needed to construct y_hat
%we also get the residuals in step 5 by subtracting y_hat from the original
%signal.
l2 = 10;
for itr = 1:totalcells
    sum = 0;
    for itr2 = 1:l2
        sum = sum + toepstore{itr,1}(1:160,itr2)*fillstore{itr,1}(itr2,1);
    end
    y_hat{itr,1} = sum;
    e2{itr,1} = C{itr,1} - y_hat{itr,1};
end
%create column vectors to store the residuals from steps 4 and 5
estep4 = cell2mat(e1);
estep5 = cell2mat(e2);
MSE_compare_4_5 = (estep4 - estep5).'*(estep4 - estep5)/length(estep4);
%%
%STEP 6 Reconstructing the Y from the residuals found previously
y_reconstructed4 = cell(totalcells,1);
y_reconstructed5 = cell(totalcells,1);
%this for loop reconstructs the y signal using the residuals calculated in
%steps 4 and 5, it shouldn't matter which we use because the residuals are
%virtually the same
for itr = 1:totalcells
    y_reconstructed4{itr,1} = y_hat{itr,1} + e1{itr,1};
    y_reconstructed5{itr,1} = y_hat{itr,1} + e2{itr,1};
end
y_recon4 = cell2mat(y_reconstructed4);
y_recon5 = cell2mat(y_reconstructed5);

MSE_recon_4 = (yremainder - y_recon4).'*(yremainder - y_recon4)/y_length;
MSE_recon_5 = (yremainder - y_recon5).'*(yremainder - y_recon5)/y_length;
%%
%STEP 7 Finding the MSE from quantizing the residuals then reconstructing
%the signal. I went back later and put in the nested for loops to gather
%multiple MSE data points for varying R and Alpha values.
k = 1;                                
j = 160;
Cresid = cell(totalcells,1);
Cresid_thresh = cell(totalcells,1);
Cresid_thresh2 = cell(totalcells,1);
%Creating the loop for storing everything into a residual cell array, and
%one for threshold
for itr=1:totalcells    
    %For loop storing audio file into 160 sample size cells
    Cresid{itr,1} = estep4(k:j); 
    Cresid_thresh{itr,1} = estep4(k:j);
    Cresid_thresh2{itr,1} = estep4(k:j);
    %increment to move to next 160 elements in array
    k = k+160;
    j = j+160;
    
end
MSEMatrixYHat = zeros(8,7);
%nested for loop used to iterate through different quantization bit rates
%as well as different alpha values.
for r=1:1:8
    alphaindex = 1;
    L = 2^r;
    for alpha = 0.5:0.5:3.5
        Cresid_thresh = Cresid;
        %For loop quantizing the residuals
        for itr=1:totalcells
            Cresid_mean = mean(Cresid{itr,1});
            Cresid_stddev = std(Cresid{itr,1});
            %setting upper and lower bounds 1 stddev away
            %For loop eliminating outliers
            for t=1:blocksize
                upperbound = Cresid_mean + alpha*Cresid_stddev;
                lowerbound = Cresid_mean - alpha*Cresid_stddev;
                if Cresid_thresh{itr,1}(t,1) > upperbound
                    Cresid_thresh{itr,1}(t,1) = upperbound;
                elseif Cresid_thresh{itr,1}(t,1) < lowerbound
                    Cresid_thresh{itr,1}(t,1) = lowerbound;
                end
            end
        end
        for i=1:totalcells
            Cresid_q = (max(Cresid_thresh{i,1})-min(Cresid_thresh{i,1}))/(L-1);
            Cresid_yq = round(Cresid_thresh{i,1}/Cresid_q)*Cresid_q;
            Cresid_thresh{i,1} = Cresid_yq;
        end
        Cresid_yqcolumn = cell2mat(Cresid_thresh);
        Yresid_recon = cell2mat(y_hat) + Cresid_yqcolumn;
        %sound(Yresid_recon);
        MSE_resid = (yremainder-Yresid_recon).'*(yremainder-Yresid_recon)/y_length;
        MSEMatrixYHat(r,alphaindex) = MSE_resid;
        alphaindex = alphaindex + 1;
    end
end
%sound(Yresid_recon);
%%
%STEP 8 Creating the different graphes needed for the plots. These are the
%3D model graphes. The 2D plots were done in the command line using the
%scatter function, such as "scatter(0.5:0.5:3.5, MSEMatrixYQ(8,:))" for
%example

%MSE_resid = (yremainder-Yresid_recon).'*(yremainder-Yresid_recon)/y_length;

% scatter([4 8 9], MSEMatrix_step_9(:,6))
% xlabel('r')
% ylabel('MSE')
% title('MSE Reconstructed from Coefficients and Residuals')
% 
% scatter(1:1:8, MSEMatrixYQ(:,6))
% xlabel('r')
% ylabel('MSE')
% title('MSE from Directly Quantized Y')
% 
% scatter(0.5:0.5:3.5, MSEMatrixYQ(8,:))
% xlabel('alpha')
% ylabel('MSE')
% title('MSE from Directly Quantized Y')
% 
% scatter(0.5:0.5:3.5, MSEMatrixYHat(8,:))
% xlabel('alpha')
% ylabel('MSE')
% title('MSE from Quantized YHat')
% 
% scatter(1:1:8, MSEMatrixYHat(:,6))
% xlabel('row')
% ylabel('MSE')
% title('MSE from Quantized Yhat')
% 
% scatter(0.5:0.5:3.5, MSEMatrix_step_9(3,:))
% xlabel('alpha')
% ylabel('MSE')
% title('MSE from Recon Y from quantized fillcoeff and residuals')
% 
% scatter([4 8 9], MSEMatrix_step_9(:,6))
% xlabel('r')
% ylabel('MSE')
% title('MSE from Recon Y from quantized fillcoeff and residuals')
% 
% scatter(1:1:8, MSEMatrixYHat(:,6))
% xlabel('r')
% ylabel('MSE')
% title('MSE from Quantized Yhat')

MSE_resid_single = MSEMatrixYHat(1,1);
surf(0.5:0.5:3.5, 1:1:8, MSEMatrixYQ)
xlabel('alpha')
ylabel('r')
zlabel('MSE')
title('MSE from directly quantized Y')

surf(0.5:0.5:3.5, 1:1:8, MSEMatrixYHat);
xlabel('alpha')
ylabel('r')
zlabel('MSE')
title('MSE from quantized YHat')

%%
%Step 9: Quantizing the filter coefficients as well as quantizing the
%residuals found in step 7, and the reconstructing the signal with a
%specific set of R values.
fillthresh = fillstore;
e_from_aq = cell(totalcells,1);
aq = cell(totalcells,1);

previous = zeros(blocksize,1);
l2 = 10;
MSEMatrix_step_9 = zeros(3,5);
%For loop constructed similarly to the rest of the nested for loops for
%finding the MSE values for varying alpha and R values.
for r = [4 8 9]
    L = 2^r;
    y_q_hat = cell(totalcells,1);
    alphaindex = 1;
    for alpha = 0.5:0.5:3.5
        fillthresh = fillstore;
        for itr = 1:totalcells
            Amean = mean(fillthresh{itr,1});
            Astddev = std(fillthresh{itr,1});
            for i = 1:length(fillthresh{itr,1})
                Aupper = Amean + alpha*Astddev;
                Alower = Amean - alpha*Astddev;
                if fillthresh{itr,1}(i,1) > Aupper
                    fillthresh{itr,1}(i,1) = Aupper;
                elseif fillthresh{itr,1}(i,1) < Alower
                    fillthresh{itr,1}(i,1) = Alower;
                end
            end
            q = (max(fillthresh{itr,1})-min(fillthresh{itr,1}))/(L-1);
            aq{itr,1} = round(fillthresh{itr,1}/q)*q;
            fillthresh{itr,1} = aq{itr,1};
        end
        for itr1 = 1:totalcells
            sum = 0;
                for itr2 = 1:l2
                    sum = sum + toepstore{itr1,1}(1:160,itr2)*fillthresh{itr1,1}(itr2,1);
                end
            y_q_hat{itr1,1} = sum;
        end
        y_resid_hat_column = cell2mat(y_q_hat);
        Cresid_thresh = Cresid;
        for itr=1:totalcells
            Cresid_mean = mean(Cresid{itr,1});
            Cresid_stddev = std(Cresid{itr,1});
            %setting upper and lower bounds 1 stddev away
            %For loop eliminating outliers
            for t=1:blocksize
                upperbound = Cresid_mean + alpha*Cresid_stddev;
                lowerbound = Cresid_mean - alpha*Cresid_stddev;
                if Cresid_thresh{itr,1}(t,1) > upperbound
                    Cresid_thresh{itr,1}(t,1) = upperbound;
                elseif Cresid_thresh{itr,1}(t,1) < lowerbound
                    Cresid_thresh{itr,1}(t,1) = lowerbound;
                end
            end
        end
        
        %Cresid_yq = cell(totalcells,1);
        for i=1:totalcells
            Cresid_q = (max(Cresid_thresh{i,1})-min(Cresid_thresh{i,1}))/(L-1);
            Cresid_yq = round(Cresid_thresh{i,1}/Cresid_q)*Cresid_q;
            Cresid_thresh{i,1} = Cresid_yq;
        end
        %After quantizing coefficients and residuals, reconstruc the Y from
        %the quantized values and then store the MSE values into a matrix
        Cresid_yqcolumn = cell2mat(Cresid_thresh);
        y_recon_q = y_resid_hat_column + Cresid_yqcolumn;
        MSE_step_9 = (yremainder-y_recon_q).'*(yremainder-y_recon_q)/y_length;
        MSEMatrix_step_9(floor(r/3),alphaindex) = MSE_step_9;
        alphaindex = alphaindex + 1;
    end
end
surf(0.5:0.5:3.5, [4 8 9], MSEMatrix_step_9);
xlabel('alpha')
ylabel('r')
zlabel('MSE')
title('MSE from recon Y from quantized fillcoeff and residuals')


