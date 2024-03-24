function runcolourdippers

% single script to run colour dipper experiment
% choose participant ID and which experiment to run
% DHB 26/8/22
% updated 25/9/22 to optimise contrast levels for each participant & condition

close all;

subj = menu('Select participant','P1','P2','P3','Test');
WaitSecs(0.5);
expt = menu('Choose experiment','Achromatic','Red/Green','Blue/Yellow');

useVSG = 1;
useRTS = 0;   % this is glitchy so best avoided

subnamelist = {'P1','P2','P3','XX'};
exptnamelist = {'AC','RG','BY'};

E.exptpath = strcat(pwd, '\'); % Same as above, JTM
KbName('UnifyKeyNames');
Screen('Preference', 'SkipSyncTests', 1);
rng('shuffle');         % reseed the random number generator from the clock

% add the colour toolbox to the path, load in the cone fundamentals and display gamut
addpath('toolbox_v1.3');
load fundamentals_ss.mat;    % stockman-sharp 2 degree cone fundamentals
load IiyamaSpectraIrradCC.mat;          % spectral readings for the Iiyama monitor in counts
phosphors = phosphors .* 1000000;
maxcont = [1 0.1 0.88];     % maximum displayable cone contrast for each stimulus type (depends on monitor gamut)
load gammavals.mat;             % load in gamma parameters from calibration

ST.npixelsperdegree = 48;
ST.gratingsizedeg = 3;
ST.gratingsizepix = ST.gratingsizedeg * ST.npixelsperdegree;
ST.gausssigma = ST.npixelsperdegree;
ST.SF = 1;
ST.ncycles = ST.SF * ST.gratingsizedeg;
ST.duration = 0.2;
ST.ISI = 0.4;
ST.postresponseduration = 0.6;
ST.phaselist = [0 90 180 270];

E.pedcontrastsdB(1,:) = -12:6:30;   
E.pedcontrastsdB(2,:) = -6:6:36;    % for the colour conditions this is % of maximum displayable cone contrast
E.pedcontrastsdB(3,:) = -6:6:36;    
E.pedcontrastsC = (10.^(E.pedcontrastsdB./20))./100;
E.pedcontrastsC(:,1) = 0;

nSCs = [4 8 8 8 8 8 8 8];

if ~exist('Results','dir')
    mkdir('Results');
end
if ~exist(strcat(E.exptpath,'Results\', subnamelist{subj}, '\'),'dir')
    mkdir(strcat(E.exptpath,'Results\', subnamelist{subj}, '\'));
end

% here generate or load in a results file to keep track of progress
if exist(strcat('Results\',subnamelist{subj},'conditionorder.mat'),'file')
    load(strcat('Results\',subnamelist{subj},'conditionorder.mat'));
    disp('Reloaded participant settings');
else
    
    thetalist = [0 160 90];      % generic values
    if exist(strcat('Results\',subnamelist{subj},'isosettings.mat'),'file')
        load(strcat('Results\',subnamelist{subj},'isosettings.mat'),'allsettings');
        thetalist(2:3) = mean(allsettings');
        disp('Loaded isoluminant settings');
    else
        disp('Using default isoluminant values');
    end
    
    % calls a function that works out the maximum displayable contrast for
    % this participant's isoluminance settings
    maxcont(2:3) = getmaxcont(thetalist,phosphors,fundamentals)
    
    condorder(1,:) = [randperm(length(E.pedcontrastsC)) randperm(length(E.pedcontrastsC)) randperm(length(E.pedcontrastsC))];
    condorder(2,:) = [randperm(length(E.pedcontrastsC)) randperm(length(E.pedcontrastsC)) randperm(length(E.pedcontrastsC))];
    condorder(3,:) = [randperm(length(E.pedcontrastsC)) randperm(length(E.pedcontrastsC)) randperm(length(E.pedcontrastsC))];
    currentrep(1:3) = 0;
    
    save(strcat('Results\',subnamelist{subj},'conditionorder.mat'),'thetalist','condorder','currentrep','maxcont');
end

% determine the condition for this block
pedlevel = condorder(expt,currentrep(expt)+1);

bg_lms = rgb2lms(phosphors,fundamentals,[0.5; 0.5; 0.5]);

theta = thetalist(expt) * pi/180;
switch expt
    case 1
        lms = [1; 1; 1];
        lmsdiff = [1; 1; 1];
    case 2
        lms = [cos(theta); sin(theta); 0];
        lmsdiff = 2.*(lms2rgb(phosphors,fundamentals,bg_lms + (maxcont(expt).*lms.*bg_lms)) - 0.5);
    case 3
        lms = [cos(theta)/sqrt(2); cos(theta)/sqrt(2); sin(theta)];
        lmsdiff = 2.*(lms2rgb(phosphors,fundamentals,bg_lms + (maxcont(expt).*lms.*bg_lms)) - 0.5);
end


try
    
    InitializePsychSound(1);
    PPA.gb = PsychPortAudio('Open',[],[],[],[],1);
    PsychPortAudio('FillBuffer', PPA.gb, MakeBeep(440*sqrt(2),0.05,44100).*0.5);
    PPA.bb = PsychPortAudio('Open',[],[],[],[],1);
    PsychPortAudio('FillBuffer', PPA.bb, MakeBeep(440/sqrt(2),0.05,44100).*0.9);
    PPA.tc = PsychPortAudio('Open',[],[],[],[],1);
    PsychPortAudio('FillBuffer', PPA.tc, MakeBeep(440,0.05,44100).*0.9);
    PPA.tclong = PsychPortAudio('Open',[],[],[],[],1);
    longbeep = [MakeBeep(440,0.05,44100).*0.9 MakeBeep(440,ST.duration+ST.ISI-0.05,44100).*0 MakeBeep(440,0.05,44100).*0.9];
    PsychPortAudio('FillBuffer', PPA.tclong, longbeep);
    
    oldVisualDebugLevel = Screen('Preference', 'VisualDebugLevel', 3);
    oldSuppressAllWarnings = Screen('Preference', 'SuppressAllWarnings', 1);
    screens=Screen('Screens');
    screenNumber=max(screens);
    
    if useVSG  % set up visage hardware
        
        crsStartup;
        CheckCard = vsgInit;
        if (CheckCard < 0)
            disp('VSG did not initialise');
            return
        end
        
        dc = crsGetSystemAttribute(CRS.DEVICECLASS);
        if dc==7
            disp('ViSaGe initialised');
        else
            disp('No ViSaGe detected');
            return;
        end
        
        crsSet42bitColourMode;
        
        width = crsGetScreenWidthPixels; 		% returns screen width
        height = crsGetScreenHeightPixels;		% returns screen height
        crsSetSpatialUnits(CRS.PIXELUNIT);				% do everything in pixels, so I know what's going on
        ST.width = width;
        ST.height = height;
        ST.framerate = round(crsGetSystemAttribute(CRS.FRAMERATE));		% in Hz
        nframestodisplay = ST.duration*ST.framerate;
        nframestowait = ST.ISI*ST.framerate;
        ST.greylevel = 0.5;
        
        crsSetDrawOrigin(width/2, height/2);	% set origin to the middle of a framestore half page
        totalpages = crsGetSystemAttribute(CRS.NUMVIDEOPAGES);
        for b = 1:totalpages
            crsClearPage(b, 0);
        end
        
        blanktexture = ones(height,width,3).*ST.greylevel;
        blanktexture = doRGBgamma(blanktexture,allp);
        crsSetDrawPage(CRS.VIDEOPAGE, 1, 0);
        crsDrawMatrix42bitColour(0, 0, blanktexture); 	% load to framestore
        crsSetDisplayPage(1);
        crsFrameSync;
        
        % also open a PTB window to blank the other screen and create the fusion lock
        rect = [1 1 1024 1024];
        [w, winRect] = Screen('OpenWindow', screenNumber, 0);  %, rect);
        Screen('FillRect', w, 0);
        Screen('Flip', w);
        
            HideCursor;
    else    % use psychtoolbox window
        
        rect = [1 1 1024 1024];
        [w, winRect] = Screen('OpenWindow', screenNumber, 128, rect);
        ST.greylevel = 0.5;
        
        Screen('FillRect', w, 255*ST.greylevel);
        Screen('Flip', w);
        
        [width, height] = Screen('WindowSize', w);
        ST.width = width;
        ST.height = height;
        ifi = Screen('GetFlipInterval', w);
        ifims = ifi * 1000 * 2; % Doubled because of frame interleaving
        
        r1 = [1 1 ST.gratingsizepix ST.gratingsizepix];
        destRectL = CenterRectOnPoint(r1, width*0.25, height*0.5);
        destRectR = CenterRectOnPoint(r1, width*0.75, height*0.5);
        
    end

    borderanglelist = randperm(360) .* pi./180;
    borderradius = ST.npixelsperdegree*3;
    ST.borderxpos = sin(borderanglelist)*borderradius;
    ST.borderypos = cos(borderanglelist)*borderradius;
    borderanglelist = randperm(360) .* pi./180;
    borderradius = ST.npixelsperdegree*3.5;
    ST.borderxpos = [ST.borderxpos sin(borderanglelist)*borderradius];
    ST.borderypos = [ST.borderypos cos(borderanglelist)*borderradius];
    borderanglelist = randperm(360) .* pi./180;
    borderradius = ST.npixelsperdegree*4;
    ST.borderxpos = [ST.borderxpos sin(borderanglelist)*borderradius];
    ST.borderypos = [ST.borderypos cos(borderanglelist)*borderradius];
    
    ST.borderlevels(1,:) = round(255.*rand(1,length(ST.borderxpos)));
    ST.borderlevels(2,:) = round(255.*rand(1,length(ST.borderxpos)));
    ST.borderlevels(3,:) = round(255.*rand(1,length(ST.borderxpos)));
    
    if useVSG
        
        % gamma correct the fusion lock squares
        borderlevsforgamma(1,:,:) = ST.borderlevels';
        borderlevsGC = doRGBgamma(borderlevsforgamma./255,allp);
        ST.borderlevels = 255.*squeeze(borderlevsGC)';
        
        % make the fusion lock image on a black background (black is transparent in the overlay)
        Screen('FillRect', w, 0);
        Screen('DrawDots', w, [(ST.width/2)+ST.borderxpos+(ST.width/4); (ST.height/2)+ST.borderypos], 8, ST.borderlevels, [], 0);
        Screen('DrawDots', w, [(ST.width/2)+ST.borderxpos-(ST.width/4); (ST.height/2)+ST.borderypos], 8, ST.borderlevels, [], 0);
        
        Screen('DrawLine', w, [1 1 1], (ST.width/4-4)-1, (ST.height/2)-1, (ST.width/4+5)-1, (ST.height/2)-1, 2);
        Screen('DrawLine', w, [1 1 1], (ST.width/4)-1, (ST.height/2+4)-1, (ST.width/4)-1, (ST.height/2-4)-1, 2);
        Screen('DrawLine', w, [1 1 1], (3*ST.width/4-4)-1, (ST.height/2)-1, (3*ST.width/4+5)-1, (ST.height/2)-1, 2);
        Screen('DrawLine', w, [1 1 1], (3*ST.width/4)-1, (ST.height/2+4)-1, (3*ST.width/4)-1, (ST.height/2-4)-1, 2);
        
        rectforframe = [0 0 ST.width ST.height];
        fusionim = Screen('GetImage',w,rectforframe,'drawBuffer');
        s = size(fusionim);
        fusionim = imresize(fusionim,[s(1) s(2)/2],'nearest');  % resize with no interpolation
        
        crsSetDrawPage(CRS.OVERLAYPAGE, 1, 0);
        crsSetPen1([1 0 0]);
        crsDrawMatrix42bitColour(0, 0, double(fusionim)./255);
        crsSetCommand(CRS.OVERLAYMASKMODE);
        
        Screen('FillRect', w, 0);
        Screen('Flip',w);
    end
    
    window = make_soft_window(ST.gratingsizepix,ST.gratingsizepix,0.9);
    for p = 1:length(ST.phaselist)
        stimulus{p} = mkgrating(ST.gratingsizepix,ST.ncycles, 90, ST.phaselist(p), 1);
    end
    
    SC = getstaircasestruct(nSCs(pedlevel),100*E.pedcontrastsC(expt,pedlevel));
    R.ntrials = zeros(SC.ncases,length(SC.levels));
    R.ncorrect = R.ntrials;
    
    if useVSG
        crsSetDisplayPage(1);
        crsFrameSync;
    else
        Screen('FillRect', w, 255*ST.greylevel);
        drawfixation(w,ST);
        lastflip = Screen('Flip', w);
    end
    
    PsychPortAudio('Start', PPA.gb);
    mouseloop;
    
    if useVSG
        crsSetDisplayPage(1);
        crsFrameSync;
    else
        Screen('FillRect', w, 255*ST.greylevel);
        drawfixation(w,ST);
        lastflip = Screen('Flip', w);
    end
    
    alltrials.subj = subnamelist{subj};
    alltrials.expt = exptnamelist{expt};
    tic;
    
    trialcounter = 0;
    while sum(SC.finished)<SC.ncases		% do trials if staircases are incomplete
        
        trialcounter = trialcounter + 1;
        testinterval = ceil(rand*2);      % test in interval 1 or 2
        trialphase = ceil(rand*length(ST.phaselist));
        
        alltrials.testinterval(trialcounter) = testinterval;
        alltrials.phase(trialcounter) = ST.phaselist(trialphase);
        
        currenttrial = 0;
        while currenttrial==0                   % select current condition
            trialcond = ceil(rand*SC.ncases);
            if SC.finished(trialcond)==0
                currenttrial = trialcond;
            end
        end
        
        alltrials.staircase(trialcounter) = currenttrial;
        alltrials.pedcontrast(trialcounter) = 100*E.pedcontrastsC(expt,pedlevel);
        alltrials.targetcontrast(trialcounter) = 100*SC.levelsM(currenttrial,SC.currentno(currenttrial));
        
        % generate stimuli
        switch currenttrial
            case 1      % monocular left
                TLstim = stimulus{trialphase}.*window.*(E.pedcontrastsC(expt,pedlevel)+SC.levelsM(currenttrial,SC.currentno(currenttrial)));      % target interval, left eye
                TRstim = stimulus{trialphase}.*0;                                           % target interval, right eye
                PLstim = stimulus{trialphase}.*window.*E.pedcontrastsC(expt,pedlevel);           % null interval, left eye
                PRstim = stimulus{trialphase}.*0;                                           % null interval, right eye
            case 2      % monocular right
                TLstim = stimulus{trialphase}.*0;                            % target interval, right eye
                TRstim = stimulus{trialphase}.*window.*(E.pedcontrastsC(expt,pedlevel)+SC.levelsM(currenttrial,SC.currentno(currenttrial)));      % target interval, left eye
                PLstim = stimulus{trialphase}.*0;                            % null interval, right eye
                PRstim = stimulus{trialphase}.*window.*E.pedcontrastsC(expt,pedlevel);      % null interval, left eye
            case {3,4}  % binocular
                TLstim = stimulus{trialphase}.*window.*(E.pedcontrastsC(expt,pedlevel)+SC.levelsM(currenttrial,SC.currentno(currenttrial)));      % target interval, left eye
                TRstim = stimulus{trialphase}.*window.*(E.pedcontrastsC(expt,pedlevel)+SC.levelsM(currenttrial,SC.currentno(currenttrial)));      % target interval, left eye
                PLstim = stimulus{trialphase}.*window.*E.pedcontrastsC(expt,pedlevel);      % null interval, left eye
                PRstim = stimulus{trialphase}.*window.*E.pedcontrastsC(expt,pedlevel);      % null interval, left eye
            case 5      % hbin left
                TLstim = stimulus{trialphase}.*window.*(E.pedcontrastsC(expt,pedlevel)+SC.levelsM(currenttrial,SC.currentno(currenttrial)));      % target interval, left eye
                TRstim = stimulus{trialphase}.*window.*E.pedcontrastsC(expt,pedlevel);                            % target interval, right eye
                PLstim = stimulus{trialphase}.*window.*E.pedcontrastsC(expt,pedlevel);                            % null interval, left eye
                PRstim = stimulus{trialphase}.*window.*E.pedcontrastsC(expt,pedlevel);                            % null interval, right eye
            case 6      % hbin right
                TLstim = stimulus{trialphase}.*window.*E.pedcontrastsC(expt,pedlevel);                            % target interval, left eye
                TRstim = stimulus{trialphase}.*window.*(E.pedcontrastsC(expt,pedlevel)+SC.levelsM(currenttrial,SC.currentno(currenttrial)));       % target interval, right eye
                PLstim = stimulus{trialphase}.*window.*E.pedcontrastsC(expt,pedlevel);                            % null interval, left eye
                PRstim = stimulus{trialphase}.*window.*E.pedcontrastsC(expt,pedlevel);                            % null interval, right eye
            case 7      % dich left
                TLstim = stimulus{trialphase}.*window.*SC.levelsM(currenttrial,SC.currentno(currenttrial));    % target interval, left eye
                TRstim = stimulus{trialphase}.*window.*E.pedcontrastsC(expt,pedlevel);                           % target interval, right eye
                PLstim = stimulus{trialphase}.*0;                                                 % null interval, left eye
                PRstim = stimulus{trialphase}.*window.*E.pedcontrastsC(expt,pedlevel);                           % null interval, right eye
            case 8      % dich right
                TLstim = stimulus{trialphase}.*window.*E.pedcontrastsC(expt,pedlevel);                            % target interval, right eye
                TRstim = stimulus{trialphase}.*window.*SC.levelsM(currenttrial,SC.currentno(currenttrial));     % target interval, left eye
                PLstim = stimulus{trialphase}.*window.*E.pedcontrastsC(expt,pedlevel);                            % null interval, right eye
                PRstim = stimulus{trialphase}.*0;                                                  % null interval, left eye
        end
        
        imgMat(:,:,1) = TLstim .* lmsdiff(1,1);
        imgMat(:,:,2) = TLstim .* lmsdiff(2,1);
        imgMat(:,:,3) = TLstim .* lmsdiff(3,1);
        TLstim = (imgMat+1)./2;
        
        imgMat(:,:,1) = TRstim .* lmsdiff(1,1);
        imgMat(:,:,2) = TRstim .* lmsdiff(2,1);
        imgMat(:,:,3) = TRstim .* lmsdiff(3,1);
        TRstim = (imgMat+1)./2;
        
        imgMat(:,:,1) = PLstim .* lmsdiff(1,1);
        imgMat(:,:,2) = PLstim .* lmsdiff(2,1);
        imgMat(:,:,3) = PLstim .* lmsdiff(3,1);
        PLstim = (imgMat+1)./2;
        
        imgMat(:,:,1) = PRstim .* lmsdiff(1,1);
        imgMat(:,:,2) = PRstim .* lmsdiff(2,1);
        imgMat(:,:,3) = PRstim .* lmsdiff(3,1);
        PRstim = (imgMat+1)./2;
        
        TLstim(find(TLstim<0)) = 0;
        TLstim(find(TLstim>1)) = 1;
        TRstim(find(TRstim<0)) = 0;
        TRstim(find(TRstim>1)) = 1;
        PLstim(find(PLstim<0)) = 0;
        PLstim(find(PLstim>1)) = 1;
        PRstim(find(PRstim<0)) = 0;
        PRstim(find(PRstim>1)) = 1;
        
        % load stimuli to graphics memory
        if useVSG        % load stimuli to framestore
            
            % resize in half vertically for later horizontal expansion
            TLstim = imresize(TLstim,[ST.gratingsizepix ST.gratingsizepix/2]);
            TRstim = imresize(TRstim,[ST.gratingsizepix ST.gratingsizepix/2]);
            PLstim = imresize(PLstim,[ST.gratingsizepix ST.gratingsizepix/2]);
            PRstim = imresize(PRstim,[ST.gratingsizepix ST.gratingsizepix/2]);
            
            % apply gamma correction
            TLstim = doRGBgamma(TLstim,allp);
            TRstim = doRGBgamma(TRstim,allp);
            PLstim = doRGBgamma(PLstim,allp);
            PRstim = doRGBgamma(PRstim,allp);
            
            crsSetDrawPage(CRS.VIDEOPAGE, 1+(2*testinterval), 0);
            crsDrawMatrix42bitColour(0, 0, blanktexture); 	% load to framestore
            crsDrawMatrix42bitColour(-width/4, 0, TLstim); 	% load to framestore
            crsDrawMatrix42bitColour(width/4, 0, TRstim); 	% load to framestore
            
            crsSetDrawPage(CRS.VIDEOPAGE, 1+(2*(3-testinterval)), 0);
            crsDrawMatrix42bitColour(0, 0, blanktexture); 	% load to framestore
            crsDrawMatrix42bitColour(-width/4, 0, PLstim); 	% load to framestore
            crsDrawMatrix42bitColour(width/4, 0, PRstim); 	% load to framestore
            
        else            % load stimuli to psychtoolbox textures
            
            if testinterval==1
                I1Ltexture = Screen('MakeTexture', w, TLstim, [], [], 2);
                I1Rtexture = Screen('MakeTexture', w, TRstim, [], [], 2);
                I2Ltexture = Screen('MakeTexture', w, PLstim, [], [], 2);
                I2Rtexture = Screen('MakeTexture', w, PRstim, [], [], 2);
            else
                I2Ltexture = Screen('MakeTexture', w, TLstim, [], [], 2);
                I2Rtexture = Screen('MakeTexture', w, TRstim, [], [], 2);
                I1Ltexture = Screen('MakeTexture', w, PLstim, [], [], 2);
                I1Rtexture = Screen('MakeTexture', w, PRstim, [], [], 2);
            end
            
        end
        if useRTS
            rtsScript = getRTSscript(nframestodisplay,nframestowait);
        end
        % wait for any remaining time before starting the next trial
        WaitSecs(ST.postresponseduration-toc);
        
        % trial sequence to display stimuli
        if useVSG        % run Visage display sequence
            
            if useRTS
                crsRTSStartStream(rtsScript, CRS.SS_IMMEDIATE);
                PsychPortAudio('Start', PPA.tclong);
                WaitSecs(2*ST.duration+ST.ISI);
            else
                PsychPortAudio('Start', PPA.tclong);
                crsSetDisplayPage(3);
                crsFrameSync;
                WaitSecs(ST.duration);
                crsSetDisplayPage(1);
                crsFrameSync;
                WaitSecs(ST.ISI);
                crsSetDisplayPage(5);
                crsFrameSync;
                WaitSecs(ST.duration);
                crsSetDisplayPage(1);
                crsFrameSync;
            end
            
        else            % run PTB display sequence
            
            % interstimulus interval
            Screen('FillRect', w, 255*ST.greylevel);
            drawfixation(w,ST);
            lastflip = Screen('Flip', w, lastflip+ST.postresponseduration);
            
            % interval 1
            Screen('FillRect', w, 255*ST.greylevel);
            Screen('DrawTexture', w, I1Ltexture, [], destRectL);
            Screen('DrawTexture', w, I1Rtexture, [], destRectR);
            drawfixation(w,ST);
            lastflip = Screen('Flip', w);
            PsychPortAudio('Start', PPA.tc);
            
            % interstimulus interval
            Screen('FillRect', w, 255*ST.greylevel);
            drawfixation(w,ST);
            lastflip = Screen('Flip', w, lastflip+ST.duration);
            
            % interval 2
            Screen('FillRect', w, 255*ST.greylevel);
            Screen('DrawTexture', w, I2Ltexture, [], destRectL);
            Screen('DrawTexture', w, I2Rtexture, [], destRectR);
            drawfixation(w,ST);
            lastflip = Screen('Flip', w, lastflip+ST.ISI);
            PsychPortAudio('Start', PPA.tc);
            
            % blank screen after stimulus offset
            Screen('FillRect', w, 255*ST.greylevel);
            drawfixation(w,ST);
            lastflip = Screen('Flip', w, lastflip+ST.duration);
            
        end
        tic;
        
        breakcode = 0;
        exitcode = 0;
        while exitcode==0
            
            [x,y,buttons] = GetMouse;
            [keyIsDown, secs, keyCode] = KbCheck;
            
            if buttons(1) || keyCode(KbName('LeftArrow'))
                resp = 1;
                exitcode = 1;
            elseif buttons(end) || keyCode(KbName('RightArrow'))
                resp = 2;
                exitcode = 1;
            elseif keyCode(KbName('Escape'))
                exitcode = 1;
                resp = 0;
                breakcode = 1;
                SC.finished(:) = 1;
            end
            
        end
        
        alltrials.responsetime(trialcounter) = toc;
        tic;
        
        
        R.ntrials(currenttrial,SC.currentno(currenttrial)) = R.ntrials(currenttrial,SC.currentno(currenttrial)) + 1;		% update results structure
        if resp==testinterval
            alltrials.iscorrect(trialcounter) = 1;
            R.ncorrect(currenttrial,SC.currentno(currenttrial)) = R.ncorrect(currenttrial,SC.currentno(currenttrial)) + 1;
            SC = dostaircase(SC, currenttrial, 1);      % correct
            PsychPortAudio('Start', PPA.gb);
        else
            alltrials.iscorrect(trialcounter) = 0;
            SC = dostaircase(SC, currenttrial, 0);      % incorrect
            PsychPortAudio('Start', PPA.bb);
        end
        
        % clean up graphics memory
        if useVSG
            if useRTS
                crsRTSDestroyStream(rtsScript);
            end
        else
            Screen('Close', I1Ltexture);
            Screen('Close', I1Rtexture);
            Screen('Close', I2Ltexture);
            Screen('Close', I2Rtexture);
            
            % blank screen to time lock to the response
            Screen('FillRect', w, 255*ST.greylevel);
            drawfixation(w,ST);
            lastflip = Screen('Flip', w);
        end
        
    end
    
    
    if ~breakcode
        currentrep(expt) = currentrep(expt) + 1;
        save(strcat('Results\',subnamelist{subj},'conditionorder.mat'),'thetalist','condorder','currentrep','maxcont');
        R.levels = SC.levels;
        save(strcat(E.exptpath,'Results\', subnamelist{subj}, '\',exptnamelist{expt},num2str(pedlevel),'rep',num2str(ceil(currentrep(expt)/length(E.pedcontrastsdB))),'.mat'), 'R','alltrials','thetalist','maxcont');
    end
    
catch
    lasterr
end

Screen('CloseAll');
ShowCursor;
clear Screen;
Screen('Close');
PsychPortAudio('Close', PPA.gb);
PsychPortAudio('Close', PPA.bb);
PsychPortAudio('Close', PPA.tc);
PsychPortAudio('Close', PPA.tclong);


% export csv text file with trial-by-trial results
fid = fopen(strcat(E.exptpath,'Results\', subnamelist{subj}, '\',exptnamelist{expt},num2str(pedlevel),'rep',num2str(ceil(currentrep(expt)/length(E.pedcontrastsdB))),'.csv'),'w');
fprintf(fid,'Subject,Experiment,Condition,PedestalContrast,TargetContrast,TargetInterval,Phase,IsCorrect,ResponseTime\n');
for s = 1:length(alltrials.testinterval)
    outputvect = [alltrials.staircase(s), alltrials.pedcontrast(s), alltrials.targetcontrast(s), alltrials.testinterval(s), alltrials.phase(s), alltrials.iscorrect(s), alltrials.responsetime(s)];
    fprintf(fid,strcat(alltrials.subj, ',', alltrials.expt, ',%2.0f,%2.3f,%2.3f,%2.0f,%2.0f,%2.0f,%2.3f\n'),outputvect);
end
fclose(fid);

if useVSG
    VSGwarmup;
end

% show figure indicating progress through experiments
figure(1);
barh(3:-1:1,currentrep);
axis([0 24 0 4]);
title(strcat(subnamelist{subj},'-Progress'));
text(11,3,'Achromatic');
text(11,2,'Red/Green');
text(11,1,'Blue/Yellow');
colormap(cool);

end
%--------------------------------------------------------------------------
function output = doRGBgamma(input,allp)

% this version for RGB stimuli, stored in an X * Y * 3 matrix
% gamma corrects the stimuli before sending to Visage
% adapted from Mark's code, DHB 29.01.08

output = zeros(size(input));
for rgbplane = 1:3
    
    i = input(:,:,rgbplane);
    k = allp(rgbplane,1);
    Lmax = allp(rgbplane,2);
    j0 = allp(rgbplane,3);
    gamma = allp(rgbplane,4);
    
    i0 = 0;
    imax = 1;                                               % Bits++ always scaled between 0 and 1
    imean = (i0+imax)/2;
    jmax = 1;
    
    Lmin = k + (Lmax-k)*(max(-j0,0)/(jmax-j0) ).^gamma;     % Eqn 2, with j set to 0, to get Lmin
    Lmin = max(Lmin,0);                                     % ensure Lmin not <0
    Lmean = (Lmin+Lmax)/2;
    L = Lmean + (Lmax-Lmean)*(i-imean)/(imax-imean);        % desired luminance values Eqn 4
    j = ((L - k)/(Lmax-k)).^(1/gamma)*(jmax - j0) + j0;     % These are the gamma-corrected lut values, j: Eqn 3
    %output = max(real(j),j0);                                     % Eqn 3 conditional
    output(:,:,rgbplane) = max(j,0);
    
end


end
%--------------------------------------------------------------------------
function drawfixation(w,ST)

% creates the fusion lock and fixation markers (non-VSG mode only)

Screen('DrawDots', w, [(ST.width/2)+ST.borderxpos+(ST.width/4); (ST.height/2)+ST.borderypos], 12, ST.borderlevels, [], 0);
Screen('DrawDots', w, [(ST.width/2)+ST.borderxpos-(ST.width/4); (ST.height/2)+ST.borderypos], 12, ST.borderlevels, [], 0);

Screen('DrawLine', w, [0], (ST.width/4-4)-1, (ST.height/2)-1, (ST.width/4+4)-1, (ST.height/2)-1, 2);
Screen('DrawLine', w, [0], (ST.width/4)-1, (ST.height/2+4)-1, (ST.width/4)-1, (ST.height/2-4)-1, 2);
Screen('DrawLine', w, [0], (3*ST.width/4-4)-1, (ST.height/2)-1, (3*ST.width/4+4)-1, (ST.height/2)-1, 2);
Screen('DrawLine', w, [0], (3*ST.width/4)-1, (ST.height/2+4)-1, (3*ST.width/4)-1, (ST.height/2-4)-1, 2);

end
%--------------------------------------------------------------------------
function mouseloop

exitcode = 0;

while exitcode==0
    [x,y,buttons] = GetMouse;
    
    if sum(buttons)>0
        exitcode = 1;
    end
end

end
%--------------------------------------------------------------------------------------------------
function imag1 = mkgrating(Regionsize, f, o, p, c)

%TSM; 26.6.03
% modified by DHB to make single component gratings, scaled from -1 to 1
% f is spatial frequency, scaled as cycles per image
% o is orientation (degrees)
% p is phase (degrees relative to centre)
% c is contrast

p = p * pi / 180;
o = o * 2 * pi / 360; % Convert from degrees to radians
f = f/Regionsize;
x0 = ((Regionsize+1) / 2);
y0 = x0;

u = f .* cos(o) * 2 * pi;
v = f .* sin(o) * 2 * pi;

imag1 = zeros(Regionsize, Regionsize);
[xx, yy] = meshgrid(1:Regionsize, 1:Regionsize);

imag1(:, :) = (c .* sin(u .* (xx-x0) + v .* (yy-y0) + p));

end
%--------------------------------------------------------------------------------------------------
function h = gausswindow(n, std)

% creates a 2d gaussian window, n*n pixels, with a sigma of std

i = repmat(1-(n/2):n/2, n, 1);
j = rot90(i);
h = exp(-((i.^2) ./ (2 .* std.^2)) - ((j.^2) ./ (2 .* std.^2)));

end
%--------------------------------------------------------------------------------------------------
function SC = getstaircasestruct(ncases,pedcontrast)

% now with separate staircase levels for different conditions
% dichoptic can get up to 100% (40dB)
% other conditions up to the closest dB value to 100-pedestal

maxval = floor(20*log10(100-pedcontrast));  % maximum dB value
minval = maxval - 57;
SC.pedlevelC = pedcontrast;
SC.ncases = ncases;                % staircase pairs
SC.stepsize = 3;		% staircase step size (dB)
SC.downrule = 3;		% decreases contrast after x correct responses
SC.uprule = 1;			% increases contrast after x incorrect responses
SC.maxtrials = 70;
SC.maxreversals = 12;   %
SC.minlevel = [minval minval minval minval minval minval -17 -17];
SC.maxlevel = [maxval maxval maxval maxval maxval maxval 40 40];
for cond = 1:SC.ncases
    SC.levels(cond,:) = SC.minlevel(cond):SC.stepsize:SC.maxlevel(cond);	   % possible test contrast values
end
SC.levelsM = 10.^(SC.levels./20)/100;                  % contrast levels in michelson
SC.startpoint(1:ncases) = 14;   % starting point in 'staircase units'
for cond = 1:SC.ncases
    SC.startlev(cond) = SC.levels(cond,SC.startpoint(cond));  % starting point in dB
end
SC.ntrials(1:SC.ncases) = 0;			% reset staircase variables for this block
SC.nreversals(1:SC.ncases) = 0;
SC.finished(1:SC.ncases) = 0;
SC.nright(1:SC.ncases) = 0;
SC.nwrong(1:SC.ncases) = 0;
SC.currentno(1:SC.ncases) = SC.startpoint;		% starting point for the staircase (staircase units)
SC.lastdir(1:SC.ncases) = 0;

end
%--------------------------------------------------------------------------
function SC = dostaircase(SC, currenttrial, response)

% this adjusts the staircase structure (SC) so that it is correct for the next trial
% current trial gives the condition (here 1-4) of the trial which has just been run
% response is either 0 (incorrect) or 1 (correct)

SC.ntrials(currenttrial) = SC.ntrials(currenttrial) + 1;		% trial counter
jump = 1;
if SC.nreversals(currenttrial)<2    % after the first trial nreverals = 1, so setting this to 2 means we really start collecting on the 1st reversal
    jump = 2;                                                   % goes in bigger steps if we're before the first reversal
    SC.ntrials(currenttrial) = SC.ntrials(currenttrial) - 1;    % don't count trials before the first reversal
end

if response==0					% add to the numbers of right or wrong responses at this level
    SC.nwrong(currenttrial) = SC.nwrong(currenttrial) + 1;
    % SC.nright(currenttrial) = 0;
else
    SC.nright(currenttrial) = SC.nright(currenttrial) + 1;
    % SC.nwrong(currenttrial) = 0;
end

thisdir = SC.lastdir(currenttrial);
if SC.nwrong(currenttrial)==SC.uprule			% need to increment
    SC.currentno(currenttrial) = SC.currentno(currenttrial) + jump;
    SC.nwrong(currenttrial) = 0;                    % reset to 0
    thisdir = 1;									% going up
elseif SC.nright(currenttrial)==SC.downrule		% need to decrement
    SC.currentno(currenttrial) = SC.currentno(currenttrial) - jump;
    SC.nright(currenttrial) = 0;                    % reset to 0
    thisdir = -1;									% going down
end

if thisdir~=SC.lastdir(currenttrial)		% have we changed direction?
    SC.nreversals(currenttrial) = SC.nreversals(currenttrial) + 1;
    SC.lastdir(currenttrial) = thisdir;
end

if SC.currentno(currenttrial)>length(SC.levels)			% check to see if we've gone too high
    SC.currentno(currenttrial) = length(SC.levels);
elseif SC.currentno(currenttrial)<1						% or too low
    SC.currentno(currenttrial) = 1;
end

if SC.ntrials(currenttrial)>=SC.maxtrials				% has the staircase finished?
    SC.finished(currenttrial) = 1;
elseif SC.nreversals(currenttrial)>=SC.maxreversals
    SC.finished(currenttrial) = 1;
end

end
%--------------------------------------------------------------------------
function imLG = makeloggabor(imSize,f0,theta0,omega,h,phi,logGabType)

% imLG = makeloggabor(imSize,f0,theta0,omega,h,logGabType)
% Produces either a Cartesian- or a polar-separable log-Gabor element.
% Input args:    imSize = width and height of output image
%    f0 = centre spatial frequency in cycles/image
%    theta0 = centre orientation in degrees
%    omega = spatial frequency bandwidth in octaves (FWHH)
%    h = orientation bandwidth in degrees (HWHH)
%    phi = phase angle in degrees
%    logGabType = `c’ or `p’ for Cartesian/polar-separable
% Output:    imLG = log Gabor filter element (imSize x imSize)
% See Baker, Summers, Baldwin & Meese (2022), PLoS ONE, doi: 10.1371/journal.pone.0267056
% Distributed under the Creative Commons Attribution-ShareAlike 4.0 International license

theta0 = theta0*pi/180; % convert all angular parameters to radians

h = h*pi/180;

phi = phi*pi/180;

u = meshgrid(1:imSize,1:imSize)-((imSize+2)/2); % set up coordinates

v = u';

f = sqrt(u.^2 + v.^2); % radial (spatial frequency) coordinate

theta = atan2(v,u); % angular (orientation) coordinate

uft = f.* cos(theta-theta0);

switch logGabType
    
    case 'c' % Cartesian-separable log Gabor
        
        numer = -(log2((f.*abs(cos(theta-theta0)))./(f0)).^2);
        
        denom = 2*(0.424*omega)^2;
        
        logGab1D = exp(numer./denom); % Implementation of Eq 2
        
        k = (log2(abs(cos(h)))/(0.424*omega))^2;
        
        eta = f0*sin(h)*sqrt(1/(log(4)-k)); % Eq 4
        
        orthFunc = exp((-(f.*sin(theta-theta0)).^2)./(2 * eta^2)); % Eq 3
        
    case 'p' % polar-separable log Gabor
        
        logGab1D = exp((-(log(f./f0)).^2)./(2*log(omega)^2)); % Eq 5
        
        sinDiff = sin(theta)*cos(theta0)-cos(theta)*sin(theta0);
        
        cosDiff = cos(theta)*cos(theta0)+sin(theta)*sin(theta0);
        
        thetaDiff = abs(atan2(sinDiff,cosDiff));
        
        orthFunc = exp((-thetaDiff.^2)./(2*h^2)); % Eq 6
        
        orthFunc = orthFunc + rot90(orthFunc,2);
        
end

logGab2D = logGab1D.* orthFunc; % combine the two filter components

cx1 = ones(imSize).* complex(0,0); % adjust the log-Gabor to be

cx2 = ones(imSize).* complex(0,0); % in the requested phase

cx1(uft>0) = complex(logGab2D(uft>0).*sin(phi),-logGab2D(uft>0).*cos(phi));

cx2(uft<0) = complex(logGab2D(uft<0).*sin(phi),logGab2D(uft<0).*cos(phi));

cxLogGab2D = cx1 + cx2;

realImage = fftshift(real(ifft2(fftshift(cxLogGab2D))));

imLG = realImage./ max(abs(realImage(:)));

% Filters are individually peak-normalised in this script. If generating
% multiple filters, you may wish to remove this line and rescale to the
% global peak of the filter bank, to balance power across all filters.

end
%--------------------------------------------------------------------------
function maxcont = getmaxcont(thetalist,phosphors,fundamentals)

% works out the maximum displayable contrast for a given angle of theta

bg_lms = rgb2lms(phosphors,fundamentals,[0.5; 0.5; 0.5]);

contrasts = 0:0.001:1;
for cond = 1:2
    
    theta = thetalist(cond+1) * pi/180;  % convert to radians
    
    switch cond
        case 1
            lms = [cos(theta); sin(theta); 0];
        case 2
            lms = [cos(theta)/sqrt(2); cos(theta)/sqrt(2); sin(theta)];
    end
    
    for c = 1:length(contrasts)
    lmsdiff = lms2rgb(phosphors,fundamentals,bg_lms + (bg_lms.*lms.*contrasts(c)));
    lmsdiff = (lmsdiff.*2)-1;   % scale from -1:1
    maxdiffP(c) = max(abs(lmsdiff));
    lmsdiff = lms2rgb(phosphors,fundamentals,bg_lms - (bg_lms.*lms.*contrasts(c)));
    lmsdiff = (lmsdiff.*2)-1;   % scale from -1:1
    maxdiffN(c) = max(abs(lmsdiff));    
    
    end
    indexP = max(find(maxdiffP<1));
    indexN = max(find(maxdiffN<1));
    maxcont(cond) = contrasts(min(indexP,indexN));
    
end

end
%--------------------------------------------------------------------------------------------------
function mask = make_soft_window(W,H,D)

% Mark's code for making a raised cosine window

% SYNTAX: mask = make_soft_window(W,H,[D])
% returns an array 'mask' that is 1 inside the circular window, shading to zero outside
% W, H are the width and height of the whole (rectangular or square) array, in pixels
% Diameter of the soft window at half-height defaults to 0.90 units
%    where 1 unit = image width or height (whichever is the smaller)
% Smoothing is by convolution with a cosine half-cycle of width 0.1 units
% Optional parameter D specifies this diameter (in relative units, range 0 -> 1)
% MAG, 27.2.04

%soft window parameters
if nargin<3, D = 0.9; end % sets default diameter to 0.9
radius = min(W*D/2,H*D/2);% radius in pixels
blur = 2*(min(W/2,H/2) - radius);  % blur half-cycle
L = blur;
X1 = [-L/2:L/2];

% 1-D blur function (applied twice, in x and y)
WinKernel = cos(X1*pi/L); % half-cycle cosine
%image coordinates - X and Y arrays
X = [1:W] - W/2;
Y = [1:H] - H/2;
xx = repmat(X,H,1);
yy = repmat(Y',1,W);

% make circular soft window
mask = single((xx.*xx + yy.*yy) < radius^2); % logical 0 outside the circle,1 inside it
mask = conv2(WinKernel,WinKernel,mask,'same'); 	% smooth the mask
mask = mask/max(max(mask));						% scale the mask 0-1
% figure(2);plot(X,mask(H/2,:),'r-',Y,mask(:,W/2));
mask = double(mask);
end
%--------------------------------------------------------------------------
function rtsScript = getRTSscript(nframestodisplay,nframestowait)

% ViSaGe RTS page cycling script for accurate timing
rtsScript = crsRTSCreateStream(0);

crsRTSAddString(rtsScript,'program(displayStimulus1);');
crsRTSAddString(rtsScript,'Long nframestoshow;');
crsRTSAddString(rtsScript,'Long nframestopause;');
crsRTSAddString(rtsScript,strcat('nframestoshow = ',num2str(nframestodisplay),';'));
crsRTSAddString(rtsScript,strcat('nframestopause = ',num2str(nframestowait),';'));

% pages are zero-indexed in RTS scripts, so subtract 1 from the page no
crsRTSAddString(rtsScript,'SET_VIDEO_WINDOW(2, 0, 0);');
crsRTSAddString(rtsScript,'SUSPEND(SUS_FRAMECOUNT,nframestoshow);');
crsRTSAddString(rtsScript,'SET_VIDEO_WINDOW(0, 0, 0);');
crsRTSAddString(rtsScript,'SUSPEND(SUS_FRAMECOUNT,nframestopause);');
crsRTSAddString(rtsScript,'SET_VIDEO_WINDOW(4, 0, 0);');
crsRTSAddString(rtsScript,'SUSPEND(SUS_FRAMECOUNT,nframestoshow);');
crsRTSAddString(rtsScript,'SET_VIDEO_WINDOW(0, 0, 0);');
crsRTSAddString(rtsScript,'SUSPEND(SUS_FRAMECOUNT,1);');

crsRTSCompileStream(rtsScript,0);

end
%--------------------------------------------------------------------------
