function eMove = step_t(rawTS,srate,winStep,winSize,thresh)
% function detects horizontal eye movements

% Time is in 1 ms steps so windows are in ms, but use samples if otherwise
winSize = round(winSize); % size of test window, ms
winStep = round(winStep); % size of step between windows, ms
% winSize = round(winSize/srate); % size of test window, samples
% winStep = round(winStep/srate); % size of step, samples

wInd = 1; eMove = zeros(1,size(rawTS,2));
while 1
    % determine portion of rawTS to test
    wEnd = wInd + winSize; 
    p = round(median([wInd wEnd])); %middle point of window
    window = wInd:wEnd;
    
    prePeak = mean(rawTS(wInd:p)); % mean amplitude prior to pointer index
    postPeak = mean(rawTS(p:wEnd)); % mean amplitude after pointer index
    
    stepAmp = abs(postPeak-prePeak); % peak to peak amplitude
    
    if stepAmp > thresh
        eMove(window) = 1; 
    end
    
    wInd = wInd + winStep;
    
    if wInd + winSize > length(rawTS)
        break
    end
end

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

