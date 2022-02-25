%% PSYCHTOOLBOX - VISUAL

% fixation cross size
crosssize= 30;

% screen colours
screenColour=[255 255 255]; %black, for white [255 255 255] I think... or [1 1 1]

% font settings
fontSizeInstr = 30;
fontSizeExp = 145;
fontType='Arial';

% Open empty window
[expWindow,rect] = Screen('OpenWindow',0,screenColour);

% Hide cursor
HideCursor;

% Getting flip rate - refresh rate of graphic display
flipInterval=Screen('GetFlipInterval',expWindow);
slack = 0.75*flipInterval; % half of the flip interval used in order not to accumulate delays, half of the refreshes of the screen should be too early and half should be too late this way

% Removes the blue screen flash and minimize extraneous warnings.
oldVisualDebugLevel = Screen('Preference', 'VisualDebugLevel', 3);

% Get coordinates of the middle of screen
[mX, mY] = RectCenter(rect);

%% Fixation Cross
% Fixation cross with black background settings
FixCr=ones(crosssize,crosssize)*255; % set background
FixCr(crosssize/2-2:crosssize/2+2,:)=0; % change colour of the rest to white
FixCr(:,crosssize/2-2:crosssize/2+2)=0;

%% Prepare textures

% Prepare fixation cross texture
fixCross = Screen('MakeTexture',expWindow,FixCr);
% White screen texture
whiteScreen=Screen('MakeTexture',expWindow,255*ones(size(rect)));

% Preload draw texture and flip function for better performance
Screen('DrawTexture', expWindow, whiteScreen);
Screen('Flip',expWindow);
WaitSecs(1);

% Preload async flip
Screen('DrawTexture', expWindow, whiteScreen);
Screen('AsyncFlipBegin', expWindow);
Screen('AsyncFlipEnd', expWindow);
WaitSecs(1);

% Settings of text
Screen('TextSize', expWindow, fontSizeInstr);
Screen('TextFont', expWindow, fontType);

% Draw instruction and wait for the space bar pressing
Screen('DrawTexture', expWindow, whiteScreen);
Screen('Flip',expWindow);

% Close screen
WaitSecs(5);
Screen('CloseAll');
Screen('Preference', 'VisualDebugLevel', oldVisualDebugLevel);
