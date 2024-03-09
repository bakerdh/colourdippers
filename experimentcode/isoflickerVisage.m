function isoflickerVisage

% experiment for finding the isoluminant point
% two separate conditions: L-M and S-(L+M)
% participant adjusts the stimulus to null perceived flicker
% we define a vector in LMS space that is the ratio of
% e.g. L and M excitation, and modify the angle theta according to the user input
% this version generates stimuli on the fly using the ViSage and stereoscope
% DHB 22/7/22

close all;
tic

subj = menu('Select participant','DHB','JTM','KJH','FGS','RJH','Test');
subnamelist = {'P1','P2','P3','P4','P5','XX'};
E.subj = subnamelist{subj};

useVSG = 1;   % either use Visage or not
ntrialspercond = 10;
ntrials = ntrialspercond*2;
stimtype = 2;

E.nConds = 2;
ST.contrast = [0.1 0.86];
ST.lumcontrast = [0 0];
ST.theta{1} = 60:0.5:179;  % set of theta values for L-M
ST.theta{2} = 45:0.5:135;  % set of theta values for S-(L+M)
ST.condlist = Shuffle([ones(1,ntrialspercond) ones(1,ntrialspercond)+1]);

nISIframes = 60;
responsekeys = {'UpArrow','DownArrow'};

E.exptPath = strcat(pwd, '\');
KbName('UnifyKeyNames');
Screen('Preference', 'SkipSyncTests', 1);

% add the colour toolbox to the path, load in the cone fundamentals and display gamut
addpath('toolbox_v1.3');
load fundamentals_ss.mat;  % stockman-sharp 2 degree
load IiyamaSpectraIrradCC.mat;
phosphors = phosphors .* 1000000;
load gammavals.mat;
%ListenChar(-1);

bg_lms = rgb2lms(phosphors,fundamentals,[0.5; 0.5; 0.5]);

if ~exist('Results','dir')
    mkdir('Results');
end
% Pixels per degree at different experimental setups
ST.nPixelsPerDegree = 48;

% Grating spans 5 degrees of visual angle
ST.gratingSize = ST.nPixelsPerDegree * 3;

ST.SF = 1;
ST.nCycles = ST.SF * ST.gratingSize / ST.nPixelsPerDegree;
ST.TFtarget = 10; % Temporal frequency waveform

WaitSecs(0.01); % Important to load in the MEX file before the expt starts
GetSecs;

% Start the 'try/catch' sequence
try
    % Hardware setup
    oldVisualDebugLevel = Screen('Preference', 'VisualDebugLevel', 3);
    oldSupressAllWarnings = Screen('Preference', 'SuppressAllWarnings', 1);
    PsychGPUControl('SetDitheringEnabled', 0);
    
    screens = Screen('Screens');
    screenNumber = max(screens);
    
    if useVSG % Using a VSG
        crsStartup;
        CheckCard = vsgInit;
        if (CheckCard < 0)
            disp('VSG did not initialise');
            return
        end
        
        dc = crsGetSystemAttribute(CRS.DEVICECLASS);
        if dc==7
            %         isvisage = 1;
            disp('ViSaGe initialised');
        else
            disp('No ViSaGe detected');
            return;
        end
        
        crsSet42bitColourMode;
        
        width = crsGetScreenWidthPixels; 		% returns screen width
        height = crsGetScreenHeightPixels;		% returns screen height
        crsSetSpatialUnits(CRS.PIXELUNIT);				% do everything in pixels, so I know what's going on
        
        greylevel = 0.5;
        blanktexture = ones(width,height,3).*greylevel;
        blanktexture = doRGBgamma(blanktexture,allp);
        crsSetDrawOrigin(width/2, height/2);	% set origin to the middle of a framestore half page
        totalpages = crsGetSystemAttribute(CRS.NUMVIDEOPAGES);
        for b = 1:totalpages
            crsClearPage(b, round(greylevel*255));
        end
        
        crsSetDrawPage(CRS.VIDEOPAGE, 1, round(greylevel*255));
        crsDrawMatrix42bitColour(0, 0, blanktexture); 	% load to framestore
        crsSetDisplayPage(1);
        crsFrameSync;
        
        framerate = crsGetSystemAttribute(CRS.FRAMERATE);		% in Hz
        ifi=1000/framerate;
        ifims = ifi * 1000;
        
        % also open a PTB window to blank the other screen and create the fusion lock
        rect = [1 1 1024 1024];
        [w, winRect] = Screen('OpenWindow', screenNumber, 0); %, rect);
        Screen('FillRect', w, 0);
        Screen('Flip', w);
        ST.width = width;
        ST.height = height;
        %[ST.width, ST.height] = Screen('WindowSize', w);
        HideCursor;
    else
        rect = [0 0 1024 1024];
        [w, winRect] = Screen('OpenWindow', screenNumber, 128, rect);
        ST.greyLevel = 0.5;
        
        Screen('FillRect', w, 255*ST.greyLevel);
        Screen('Flip', w);
        
        [width, height] = Screen('WindowSize', w);
        ifi = Screen('GetFlipInterval', w);
        ifims = ifi * 1000;
        
        r1 = [1 1 ST.gratingSize ST.gratingSize];
        destRectL = CenterRectOnPoint(r1, width*0.25, height*0.5);
        destRectR = CenterRectOnPoint(r1, width*0.75, height*0.5);
        
    end
    
    ST.nFrames = round(1000 / ifims);
    
    targetWaveform = [1 -1];
    
    % Replace center of window with soft window
    window = make_soft_window(ST.gratingSize, ST.gratingSize, 0.9);
    
    % Make RGB
    window = cat(3, window, window, window);
    
    if stimtype==1
        colgrating = mkgrating(ST.gratingSize, ST.nCycles, 0, 90, 1).*0 + 1;
    elseif stimtype==2
         colgrating = mkgrating(ST.gratingSize, ST.nCycles, 90, 90, 1);
%         colgrating = makeloggabor(ST.gratingSize, ST.nCycles, 90, 0.8, 15, 90, 'c');
    else
        colgrating = circgrating(ST.gratingSize, ST.nCycles, 0, 1);
    end
    
    if ~useVSG
        nullImgMat = zeros(ST.gratingSize,ST.gratingSize,3);
        nullImgMat = (1+nullImgMat) ./ 2;
        nullTexture = Screen('MakeTexture', w, nullImgMat, [], [], 2);
    end
    
    borderanglelist = randperm(360) .* pi./180;
    borderradius = ST.nPixelsPerDegree*3;
    ST.borderxpos = sin(borderanglelist)*borderradius;
    ST.borderypos = cos(borderanglelist)*borderradius;
    borderanglelist = randperm(360) .* pi./180;
    borderradius = ST.nPixelsPerDegree*3.5;
    ST.borderxpos = [ST.borderxpos sin(borderanglelist)*borderradius];
    ST.borderypos = [ST.borderypos cos(borderanglelist)*borderradius];
    borderanglelist = randperm(360) .* pi./180;
    borderradius = ST.nPixelsPerDegree*4;
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
    
    
    % ready to start experiment, wait for some input
    
    exitcode = 0;
    while ~exitcode
        [keyIsDown, secs, keyCode] = KbCheck;
        if keyIsDown
            exitcode = 1;
        end
        
        [x, y, buttons] = GetMouse;
        if sum(buttons) > 0
            exitcode = 1;
        end
    end
WaitSecs(1);
    % just do a flip so we have the start time recorded
    if useVSG
        crsSetDisplayPage(1);
        crsFrameSync;
    else
        Screen('FillRect', w, 255*ST.greyLevel);
        Screen('DrawTexture', w, nullTexture, [], destRectL);
        % Screen('DrawDots', w, ST.dotCoords, ST.dotSize, ST.dotLevels, ST.centre, 0);  % Change 1 to 0 to draw square dots
        lastFlip = Screen('Flip', w);
    end
    
    ST.nFramesperphase = 5;
    
    allTriggerTimes = [];
    
    condtrialcount = [0 0];
    breakCode = 0;
    currentBlock = 0;
    
    trial = 0;
    while trial < ntrials
        
        trial = trial + 1;
        cond = ST.condlist(trial);
        condtrialcount(cond) = condtrialcount(cond) + 1;
        
        thetalist = ST.theta{cond};
        currentlevel = ceil(rand*length(thetalist));
        
        %startTime = lastFlip;
        phaseval = -1;
        n = 0;
        exitcode = 0;
        while ~exitcode
            
            phaseval = phaseval * -1;   % flip the stimulus phase every 5 frames
            SetMouse(width/2,round(height*(length(thetalist)-(currentlevel-1))/length(thetalist)),w);     % set the mouse pointer to the appropriate location for the random size
            
            theta = thetalist(currentlevel) * pi/180;  % convert to radians
            
            switch cond
                case 1
                    lms = [cos(theta); sin(theta); 0];
                case 2
                    lms = [cos(theta)/sqrt(2); cos(theta)/sqrt(2); sin(theta)];
            end
            
            lmsdiff = 2.*(lms2rgb(phosphors,fundamentals,bg_lms + (bg_lms.*ST.contrast(cond).*lms)) - 0.5);
            
            imgMat(:,:,1) = (colgrating .* lmsdiff(1,1) .* phaseval);
            imgMat(:,:,2) = (colgrating .* lmsdiff(2,1) .* phaseval);
            imgMat(:,:,3) = (colgrating .* lmsdiff(3,1) .* phaseval);
            
            imgMat = imgMat.*window;
            imgMat = (1+imgMat)/2;
            
            if sum(find(imgMat<0))>0
                disp('Stimulus outside gamut');
            end
            if sum(find(imgMat>1))>0
                disp('Stimulus outside gamut');
            end
            
            imgMat(find(imgMat<0)) = 0;
            imgMat(find(imgMat>1)) = 1;
            
            if useVSG
                
                imgMat = imresize(imgMat,[ST.gratingSize ST.gratingSize/2]);
                imgMat = doRGBgamma(imgMat,allp);
                
                crsSetDrawPage(CRS.VIDEOPAGE, 4+phaseval, 0);
                crsDrawMatrix42bitColour(0, 0, blanktexture); 	% load to framestore
                crsDrawMatrix42bitColour(-width/4, 0, imgMat); 	% load to framestore
                crsDrawMatrix42bitColour(width/4, 0, imgMat); 	% load to framestore
                crsSetDisplayPage(4+phaseval);
                crsFrameSync;
                clear imgMat;
            else
                
                % Add other textures here targetTexture(1,n)=, clear screen
                targetTexture = Screen('MakeTexture', w, imgMat, [], [], 2);
                
                % Function to avoid duplication below?
                Screen('FillRect', w, 255*ST.greyLevel);
                Screen('DrawTexture', w, targetTexture, [], destRectL);
                Screen('DrawTexture', w, targetTexture, [], destRectR);
                
                % Screen('DrawDots', w, ST.dotCoords, ST.dotSize, ST.dotLevels, ST.centre, 0);  % Change 1 to 0 to draw square dots
                Screen('DrawingFinished', w);
                lastFlip = Screen('Flip', w, lastFlip + ifi * (ST.nFramesperphase-0.5));
                % Screen('Close',targetTexture);
            end
            
            [keyIsDown, secs, keyCode] = KbCheck;
            if keyCode(KbName(responsekeys{1}))  % responded up
                currentlevel = currentlevel + 1;
                %                     disp(strcat('Up; level increasing to', num2str(currentlevel')));
            end
            if keyCode(KbName(responsekeys{2}))  % responded down
                currentlevel = currentlevel - 1;
                %                     disp(strcat('Down; level decreasing to', num2str(currentlevel')));
            end
            if currentlevel<1
                currentlevel = 1;
            end
            if currentlevel>length(thetalist)
                currentlevel = length(thetalist);
            end
            
            [x, y, buttons] = GetMouse;
            if sum(buttons) > 0
                exitcode = 1;
            end
            if ~keyIsDown
                currentlevel = (length(thetalist)+1)-round(length(thetalist)*(y/height));
                currentlevel = min(currentlevel, length(thetalist));
                currentlevel = max(currentlevel, 1);
            end
            
            if keyCode(KbName('.>'))  % responded quit
                exitcode = 1;
                breakCode = 1;
                trial = 10000;
            end
        end
        
        
        
        
        if useVSG
            crsSetDisplayPage(1);
            crsFrameSync;
            WaitSecs(1);
        else
            trialoffset = lastFlip;
            for frame = 1:nISIframes
                % Shows null stimuli to avoid flicker on initial frame of first trial
                
                Screen('FillRect', w, 255*ST.greyLevel);
                Screen('DrawTexture', w, nullTexture, [], destRectL);
                % Screen('DrawDots', w, ST.dotCoords, ST.dotSize, ST.dotLevels, ST.centre, 0);  % Change 1 to 0 to draw square dots
                Screen('DrawingFinished', w);
                lastFlip = Screen('Flip', w, lastFlip + ifi * 0.5);
                
            end
        end
        
        allsettings(cond,condtrialcount(cond)) = thetalist(currentlevel);
        
    end
    
    
    
    % If main try loop fails, get the last error
catch
    
    lasterr
    
end

% R.ncorrect'
% R.ntrials'
% save(strcat('Results/',E.subj,'isosettingsMRI.mat'),'R');

Screen('Preference', 'VisualDebugLevel', oldVisualDebugLevel);
Screen('Preference', 'SuppressAllWarnings', oldSupressAllWarnings);

Screen('Flip', w);
if ~useVSG
    Screen('Close', nullTexture);
end
ShowCursor;

Screen('LoadNormalizedGammaTable', w, linspace(0,1,256)' * ones(1,3), 0);

Screen('CloseAll');

if useVSG
    VSGwarmup;
end

ListenChar(1);

allsettings
mean(allsettings')

save(strcat('Results\',E.subj,'isosettings.mat'),'ST','allsettings');


thetalist = (0:1:360).*pi/180;
scale = [2 4];

for n = 1:length(thetalist)
    theta = thetalist(n);
    stimLMS = [cos(theta); sin(theta); 0]; % should be [1 -1 0]
    rgb(n,:) = lms2rgb(phosphors,fundamentals,scale(1).*stimLMS) + 0.5;
end

rgb(find(rgb<0)) = 0;
rgb(find(rgb>1)) = 1;

figure(1);
set(gcf,'Position',[-13 316 1027 481]);
subplot(1,2,1)
axis([-1 1 -1 1]);
xlabel('L');
ylabel('M');
axis square;
hold on;
plot([-1 1],[0 0],'k--');
plot([0 0],[-1 1],'k--');
for n = 1:length(thetalist)
    
    plot(cos(thetalist(n)),sin(thetalist(n)),'.','Color',rgb(n,:),'MarkerSize',20)
    
end

thetalist2 = ST.theta{1}.*pi/180;
plot(cos(thetalist2).*0.5,sin(thetalist2).*0.5,'k-','LineWidth',2);
plot(-cos(thetalist2).*0.5,-sin(thetalist2).*0.5,'k-','LineWidth',2);

for n = 1:ntrialspercond
    a = allsettings(1,n)*pi/180;
    plot(cos(a).*[-1 1],sin(a).*[-1 1],'k:','LineWidth',1);
end
a = mean(allsettings(1,:))*pi/180;
plot(cos(a).*[-1 1],sin(a).*[-1 1],'k-','LineWidth',3);

text(0.75,-0.9,E.subj);


for n = 1:length(thetalist)
    theta = thetalist(n);
    stimLMS = [cos(theta); cos(theta); sin(theta)]; % should be [1 -1 0]
    rgb(n,:) = lms2rgb(phosphors,fundamentals,scale(2).*stimLMS) + 0.5;
end

rgb(find(rgb<0)) = 0;
rgb(find(rgb>1)) = 1;

subplot(1,2,2)
axis([-1 1 -1 1]);
xlabel('L+M');
ylabel('S');
axis square;
hold on;
plot([-1 1],[0 0],'k--');
plot([0 0],[-1 1],'k--');
for n = 1:length(thetalist)
    
    plot(cos(thetalist(n)),sin(thetalist(n)),'.','Color',rgb(n,:),'MarkerSize',20)
    
end

thetalist2 = ST.theta{2}.*pi/180;
plot(cos(thetalist2).*0.5,sin(thetalist2).*0.5,'k-','LineWidth',2);
plot(-cos(thetalist2).*0.5,-sin(thetalist2).*0.5,'k-','LineWidth',2);

for n = 1:ntrialspercond
    a = allsettings(2,n)*pi/180;
    plot(cos(a).*[-1 1],sin(a).*[-1 1],'k:','LineWidth',1);
end
a = mean(allsettings(2,:))*pi/180;
plot(cos(a).*[-1 1],sin(a).*[-1 1],'k-','LineWidth',3);

text(0.75,-0.9,E.subj);

saveas(gcf,strcat('Results\',E.subj,'isosettings.jpg'),'jpeg');


toc
end
%--------------------------------------------------------------------------
function mouseloop

exitCode = 0;

while exitCode==0
    [x, y, buttons] = GetMouse;
    
    if sum(buttons) > 0
        exitCode = 1;
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
%--------------------------------------------------------------------------
function kbloop

exitCode = 0;

while exitCode==0
    [keyIsDown, secs, keyCode] = KbCheck; % Also monitor keyboard for breaks
    
    if sum(keyCode(79:80)) > 0 % Respond to left and right arrows
        exitCode = 1;
    end
end

end
%--------------------------------------------------------------------------------------------------
function mask = make_soft_window(W,H,D)

% Mark's code for making a raised cosine window

% SYNTAX: mask = make_soft_window(W,H,[D])
% ends an array 'mask' that is 1 inside the circular window, shading to zero outside
% W, H are the width and height of the whole (rectangular or square) array, in pixels
% Diameter of the soft window at half-height defaults to 0.90 units
%    where 1 unit = image width or height (whichever is the smaller)
% Smoothing is by convolution with a cosine half-cycle of width 0.1 units
% Optional parameter D specifies this diameter (in relative units, range 0 -> 1)
% MAG, 27.2.04

% Soft window parameters
if nargin < 3, D = 0.9; end % Sets default diameter to 0.9
radius = min(W*D/2, H*D/2); % Radius in pixels
blur = 2 * (min(W/2, H/2) - radius); % Blur half-cycle
L = blur;
X1 = [-L/2:L/2];

% 1-D blur function (applied twice, in x and y)
WinKernel = cos(X1*pi/L); % half-cycle cosine
% Image coordinates - X and Y arrays
X = [1:W] - W/2;
Y = [1:H] - H/2;
xx = repmat(X,H,1);
yy = repmat(Y',1,W);

% Make circular soft window
mask = single((xx.*xx + yy.*yy) < radius^2); % logical 0 outside the circle,1 inside it
mask = conv2(WinKernel,WinKernel,mask,'same'); 	% smooth the mask
mask = mask / max(max(mask));						% scale the mask 0-1
% figure(2);plot(X,mask(H/2,:),'r-',Y,mask(:,W/2));
mask = double(mask);
end
%--------------------------------------------------------------------------
function grating = circgrating(imsize, circfreq, phase, contrast)

% Function to generate circular grating.
%
% imsize is the diameter of the grating pattern in pixels
% circfreq is the frequency of the circles in Hz
% phase is the phase of the circles
% contrast is the contrast in ...

[u, v] = meshgrid((1:imsize)-(imsize+2)/2, (1:imsize)-(imsize+2)/2);
f = sqrt(u.^2 + v.^2); % Radial coordinate
theta = atan2(v, u); % Angular coordinate
f = f ./ max(f(imsize/2,:));
theta = ((theta./pi)+1)./2;

grating = cos(circfreq .* f .* 2 .* pi + phase) .* contrast;


end
%--------------------------------------------------------------------------
function Clut = makeLutPro(c,p,Nbits,ibits)

% my version of this for use with colour monitors
% DHB 7/08/07

% Syntax: Clut = makeLutPro(c,p,Nbits,ibits)
% This function creates a 3-column LUT that attenuates contrast by factor c, and applies gamma correction.
% ibits is the no of greylevel indices; Nbits is the no. of bits in the DAC
% p is the vector of display system parameters [j0, Lmax, gamma, k]
% gamma is the monitor gamma, from calibration data
% j0,Lmax, k are the other calibrated parameters of the CRT model - see document 'CRT model8.pdf'
% k is the (theoretical) luminance emitted when j=j0. Lmin (below) is the minimum luminance (at j=0)

i0 = 0;
imax = 2^ibits-1;
imean = (i0+imax)/2;
jmax = 2^Nbits-1;

for n = 1:3
    k = p(1,n) ;
    Lmax=p(2,n);
    j0=p(3,n);
    gamma=p(4,n); % from lut parameter vector, p
    
    Lmin = k + (Lmax-k)*( max(-j0,0)/(jmax-j0) ).^gamma; % Eqn 2, with j set to 0, to get Lmin
    Lmin = max(Lmin,0); % ensure Lmin not <0
    Lmean = (Lmin+Lmax)/2 ;
    i = [ i0:imax ]; % vector of greyscale values, i
    L = Lmean + c*(Lmax-Lmean)*(i-imean)/(imax-imean); % desired luminance values, scaled by contrast factor, c: Eqn 4
    j = ( (L - k)/(Lmax-k) ).^(1/gamma)*(jmax - j0) + j0; % These are the gamma-corrected lut values, j: Eqn 3
    lutval = round(max(j,j0)); % Eqn 3 conditional
    Clut(i+1,n)= lutval';   %red, green, blue columns
    Clut(1,n) = 0; % black 'reserved' colour for fixation point etc.
end

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
%--------------------------------------------------------------------------------------------------