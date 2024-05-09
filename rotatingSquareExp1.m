% rotating square experiment

% squares are made up of Gabor elements aligned with the contours

% 2 IFC MOCS expt

% reference square is actually rotating

% comparison square has stationary Gabors that are changing in phase

% task: which interval contains the faster square?

% manipulations: presentation duration and rotation speed

%

% Last Edit: GE 4/7





Screen('Preference', 'SkipSyncTests', 2);

clear; 

AssertOpenGL;

KbName('UnifyKeyNames');



debugMode = 1;

EXITNOW = 0;





try

%% SCREEN SETUP

screenid = max(Screen('Screens'));

PsychImaging('PrepareConfiguration');

PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');

[win screenRect] = PsychImaging('OpenWindow', screenid, 128);

Screen('BlendFunction', win, GL_ONE, GL_ONE);

[cx cy] = RectCenter(screenRect);

ifi = Screen('GetFlipInterval',win);

hz = 60;



white = WhiteIndex(win);

black = BlackIndex(win);

red = [255 0 0];

green = [0 255 0];

blue = [0 0 255];

gray = round(white./2);



%% EXPERIMENTAL DESIGN



nTestSpeeds = 8;

testSpeeds = linspace(1,8,nTestSpeeds)/10;  % rotational velocity (degrees/frame)

refSpeed = mean(testSpeeds);                % reference rotational velocity 

nPresTimes = 5;

presTimes = hz:hz:hz*nPresTimes;

presTimes = presTimes./1000;                % convert to ms

nRep = 20;                                  % number of repetitions of each trial



trialMat = fullfact([nPresTimes,nTestSpeeds]);  % all condition combinations

trialMat = repmat(trialMat,nRep,1);             % repeat for number of trials

nTrials = size(trialMat,1);                     % total number of trials

trialMat = trialMat(randperm(nTrials),:);       % randomize trial order



%%  STIMULUS PROPERTIES



% Gabors

si = 32;

tw = 2*si+1;

th = 2*si+1;

phase = 0;

sc = 10.0;

freq = .05;

tilt = 0;

contrast = 10.0;

aspectratio = 1.0;

gabortex = CreateProceduralGabor(win, tw, th);

texRect = [0 0 tw th];



% define gabor positions along contour of square

nGaborPerSide = 7;

nGabors = nGaborPerSide*4;

interGaborDist = 60;

xpos = linspace(cx-((nGaborPerSide-1)/2)*interGaborDist,cx+((nGaborPerSide-1)/2)*interGaborDist,nGaborPerSide);

ypos = linspace(cy-((nGaborPerSide-1)/2)*interGaborDist,cy+((nGaborPerSide-1)/2)*interGaborDist,nGaborPerSide);

gaborPos = [xpos ones(1,nGaborPerSide)*(xpos(end)+interGaborDist/2) fliplr(xpos) ones(1,nGaborPerSide)*(xpos(1)-interGaborDist/2); ...

    ones(1,nGaborPerSide)*(ypos(1)-interGaborDist/2) ypos ones(1,nGaborPerSide)*(ypos(end)+interGaborDist/2) fliplr(ypos)];

dstRects = zeros(4, nGabors);

for ii=1:nGabors

    dstRects(:, ii) = CenterRectOnPoint(texRect, gaborPos(1,ii),gaborPos(2,ii))';

end

mypars = repmat([180-phase, freq, sc, contrast, aspectratio, 0, 0, 0]',1,nGabors);

drawAngles = [90*ones(1,nGaborPerSide), zeros(1,nGaborPerSide),90*ones(1,nGaborPerSide), zeros(1,nGaborPerSide)];



% define velocity vectors for each point along square as it rotates

% do this for each of the rotational speeds that will be used

normFlow = zeros(nTestSpeeds,nGabors); % matrix storing the flow vectors

for s = 1:nTestSpeeds

    rotVel = testSpeeds(s); % degrees per frame

    rotMat = [cosd(rotVel) -sind(rotVel) ; sind(rotVel) cosd(rotVel)];

    zeroGaborPos = gaborPos;

    zeroGaborPos(1,:) = zeroGaborPos(1,:) - cx;

    zeroGaborPos(2,:) = zeroGaborPos(2,:) - cy;

    rotGaborPos = rotMat'*zeroGaborPos;

    eucDist = sqrt((zeroGaborPos(1,:)-rotGaborPos(1,:)).^2+(zeroGaborPos(2,:)-rotGaborPos(2,:)).^2);

    % eucDist is how much each element travels -- per frame -> rotational velocity

    % compute normal flow (y-component of eucDist)

    normFlow(s,:) = eucDist.*sind(rotVel);

end % describing flow for each of the test speeds



% manually adjust sign

if mod(nGaborPerSide,2)==0

    flowSign = [-ones(1,nGaborPerSide/2) ones(1,nGaborPerSide/2) ...

        -ones(1,nGaborPerSide/2) ones(1,nGaborPerSide/2) ...

        ones(1,nGaborPerSide/2) -ones(1,nGaborPerSide/2) ...

        ones(1,nGaborPerSide/2) -ones(1,nGaborPerSide/2)];

else

    flowSign = [-ones(1,(nGaborPerSide-1)/2) 0 ones(1,(nGaborPerSide-1)/2) ...

        -ones(1,(nGaborPerSide-1)/2) 0 ones(1,(nGaborPerSide-1)/2) ...

        ones(1,(nGaborPerSide-1)/2) 0 -ones(1,(nGaborPerSide-1)/2) ...

        ones(1,(nGaborPerSide-1)/2) 0 -ones(1,(nGaborPerSide-1)/2)];

end

normFlow = normFlow.*repmat(flowSign,[nTestSpeeds,1]);



% convert from pixels to phase

normFlow = normFlow.*freq*360* tw;%?



% how many Gabors are we drawing?

% everything

allGabors = 1:nGaborPerSide*4;

% only top and bottom (open figure)

topbottom = [1:7, 15:21];

% only corners

corners = [1:2, 6:9, 13:16, 20:23, 27:28];

% only sides

sides = setdiff(allGabors,corners);

% drawing condition?

drawOnly = allGabors;



% remove Gabors that we aren't drawing

dstRects = dstRects(:,drawOnly);

drawAngles = drawAngles(drawOnly);

mypars = mypars(:,drawOnly);

normFlow = normFlow(:,drawOnly);



origRects = dstRects;

origPos = gaborPos;

origAngs = drawAngles;

origPars = mypars;



ISI = 0.5;  % blank interval between stimuli (s)



%% SETUP DATA STRUCTURE

% record some information about the display and experiment in results

% structure res

res.display = Screen('Resolution',screenid);

res.trialMat = trialMat;

res.refSpeed = refSpeed;

res.testSpeeds = testSpeeds;

res.presTimes = presTimes;



% create variables for storing responses

res.refOrder = zeros(nTrials,1);            % interval order of reference and comparison

res.fasterOrder = zeros(nTrials,1);         % which one is faster (ref or comparison)

res.response = zeros(nTrials,1);            % subject's response for which is faster



%% INSTRUCTIONS

if ~debugMode % skip instructions if debugging

end



%% EXPERIMENT

if debugMode % if debugging, only show the first 5 trials

    nTrials = 10;

end



KbReleaseWait;

for t = 1:nTrials

    % determine settings for this particular trial

    this_presTime = presTimes(trialMat(t,1)); % presentation time (same for both intervals)

    this_speed = testSpeeds(trialMat(t,2)); % rotation speed

    

    % determine whether first interval contains test or reference speed

    if rand <0.5

        intSpeed(1) = this_speed;

        intSpeed(2) = refSpeed;

        res.refOrder(t) = 2;

    else

        intSpeed(1) = refSpeed;

        intSpeed(2) = this_speed;

        res.refOrder(t) = 1;

    end

    

    % which interval contains the fastest speed?

    if intSpeed(1) > intSpeed(2)

        res.fasterOrder(t) = 1;

    else

        res.fasterOrder(t) = 2;

    end

    

    % randomize rotation direction

    rotDir = logical(round(rand)); % 0 - counter-clockwise   1 - clockwise

        

    % display two intervals

    for interval = 1:2

        % determine the rotation rate of the square 

        if res.refOrder(t)==interval % if reference, spinning

            if rotDir

                rotMat = [cosd(this_speed) sind(this_speed) ; -sind(this_speed) cosd(this_speed)];

                angUpdate = -this_speed;

            else

                rotMat = [cosd(this_speed) -sind(this_speed) ; sind(this_speed) cosd(this_speed)];

                angUpdate = this_speed;

            end

            this_flow = 0;

        else % if comparison, stationary

            rotMat = [1 0 ; 0 1];

            angUpdate = 0;

            this_flow = normFlow(trialMat(t,2),:);  % flow vectors of Gabors

            if ~rotDir

                this_flow = -this_flow;

            end

        end

        % initiate starting position of Gabors

        dstRects = origRects;

        gaborPos = origPos;

        drawAngles = origAngs;

        mypars = origPars;

        

        % initiate timer

        startTime = GetSecs; 

        

        % present the stimulus

        while GetSecs-startTime < this_presTime

            Screen('DrawTextures', win, gabortex, [], dstRects, drawAngles, [], [], [], [], kPsychDontDoRotation, mypars);

            Screen('Flip',win);

        

            % update phase

            mypars(1,:) = mypars(1,:) + this_flow;

            

            % update orientation of Gabors (if actually rotating)

            drawAngles = drawAngles + angUpdate;

            

            % rotate the positions of all of the Gabors by rotMat

            gaborPos(1,:) = gaborPos(1,:) - cx;

            gaborPos(2,:) = gaborPos(2,:) - cy;

            gaborPos = rotMat'*gaborPos;

            gaborPos(1,:) = gaborPos(1,:) + cx;

            gaborPos(2,:) = gaborPos(2,:) + cy;

            for ii=1:nGabors

                dstRects(:, ii) = CenterRectOnPoint(texRect, gaborPos(1,ii),gaborPos(2,ii))';

            end  

        end % animation

        

        % ISI

        Screen('FillRect',win,gray);

        Screen('Flip',win);

        WaitSecs(ISI);

        

    end % intervals

    

    % collect response

    KbReleaseWait;

    DrawFormattedText(win,'Was the square spinning faster in the first (A) or second (L) movie?','center','center',black);

    Screen('Flip',win);

    noGoodKey = 1;

    while noGoodKey

        KbReleaseWait;

        [secs keyCode] = KbWait;

        keyPressed = KbName(keyCode);

        if strcmpi(keyPressed,'ESCAPE') % exit out of experiment

            EXITNOW = 1;

            break;

        end

        % check for valid response key

        if strcmpi(keyPressed,'a')

            res.response(t) = 1;

            noGoodKey = 0;

        elseif strcmpi(keyPressed,'l')

            res.response(t) = 2;   

            noGoodKey = 0;

        end

    end

    

    % check to see if exit key has been pressed

    if EXITNOW

        break;

    end

end

%% CLEANUP

sca;



catch

    sca;

    psychrethrow(psychlasterror);

end