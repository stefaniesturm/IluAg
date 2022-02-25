function IoA_ATrain(Subj,Block,Repetition,counterbalancing)

% Illusion of Agency experiment, Bachelor thesis of Marta Blasco and
% Master thesis of Ignacio Ruiz
% financed by TEMPLATES project (PSI2014- 52573-P, Proyectos de I+D
% “EXCELENCIA”, MINECO); RYC-2013-12577 (MINECO)
%
% TASK - Training block for acquisition trials
%
% % In the acquisition trials, the participant need to learn a mapping between two
% buttons and two sounds. They press one of two buttons after seeing a cue
% in every trial and a sound contingent to that button press is generated.

% INPUT:
% Subj, subject number, e.g. 1
% Block, block number, e.g. 1
% Repetition, 1st, 2nd, 3rd ..etc repetition of the same block, e.g., 1
% Counterbalancing: 1 or 2.

%
% RESPONSE DEVICE AND BUTTONS
% Button presses to "produce" the sounds are collected with a keyboard with two buttons (left uppercase key and left control key) pressing with two fingers of the left hand.
%
% TRIALS PER BLOCK
% There will be 20 trials per block

%
% OUTPUT: LOGFILE
% Generates a logfile for each block named “subject number_A_block number_repetition number.txt” in
% a "Res" folder in the experiment folder.
% If a block is repeated, a new logfile is created for every repeated block
% (coded in repetition number; if the same repetition number is given, it can be overwritten).
% The logfile contains the following information:
% Column 1: Trial number
% Column 2: Subject number
% Column 3: Trial Type (1 = acquisition, 2 = test)
% Column 4: Block number (1-16)
% Column 5: Repetition number
% Column 6: Cue time (response window starts)
% Column 7: Cue turns off (resp window end)
% Column 8: Button press time
% Column 9:Button press time
% Column 10: Presented sound (1 = low, 2 = high)
% Column 11: Sound time
% Column 12: cue to press RT
% Column 13: Question start time
% Column 14: Agency judgement (1-4)
% Column 15: Agency judgement time
% Column 16: Question to JoA RT
% Column 17: Condition (1 = congruent, 0 = incongruent, 2 = no press)
% Column 18: Sound-Button order (-1 = sound first; 1 = button first)
% Column 19: time of sound relative to time of press (-100 to +60)
% Column 20: Time bin we aimed at (1-8)
% Column 21: time bin we hit (1-8)

%
% SOUND STIMULI:
% 2 sounds of 100ms duration
%
% ADDITIONAL SCRIPTS AND FILES NEEDED:
% - The soundfiles in .wav format need to be placed in experiment_path\Sounds, named Ta_high_243Hz.wav and Ta_low_105Hz.wav
%


%% Set up
% --------------------------------------------------------------------
%               Paths and other general variables
%----------------------------------------------------------------------
% Paths
experiment_path = pwd;
Sounds_path = [experiment_path '\Sounds'];
Results_path = [experiment_path '\Res'];
Phase = 'ATrain' %acquisition training (ATrain), main training (MTrain), experimental block (Exp) or motor-only block (MO)

% Timing vars
RespWinDur = 0.8; % From fixation, we check for presses for 800 ms. Afterwards fixation dissapears
RTlimits = [0.41 0.59]; % RTs considered "correct" between 420 and 580 ms. If pressed earlier cross turns orange, if later (or no resp) cross turn red
Fix2QTime = [1.7];% From fixation to question (Ttrial) or to end of trial (Atrial)
ITI = [0.3:0.1:0.8];% From 0 to Fixation (Ttrial), or from end of trial to Fixation (Atrial)

%----------------------------------------------------------------------
%                       Prepare Logfile
%----------------------------------------------------------------------
mkdir(Results_path);
% creates the pathname for the logfile inside the experiment folder, e.g. logfilename = 'results\01_01_01.txt'
logfilename = [Results_path '\' sprintf('%02d',Subj) '_' sprintf('%s',Phase) '_' sprintf('%02d',Block) '_' sprintf('%02d',Repetition) '.txt'];
exists = fopen(logfilename, 'r'); % checks if logfile already exists, returns -1 if does not exist
if exists ~= -1
    overwrite = input('WARNING: Logfile already exists, overwrite? 1 = yes, 2 = no');
    if overwrite == 2
        return;
    end
end
LOGFILE = fopen(logfilename, 'w+'); % creates the logfile if nonexistent, opens it for writing and specifies
% old info will be overwritten if already exists


%--------------------------------------------------------------------
%                   Screen settings
%----------------------------------------------------------------------
backcolor = [255 255 255]; % white
res = [600 600]; % screen resolution for display window for debugging (comment for fullscreen)
textcolor = [0 0 0]; % black
textcolorcorrect = [0 255 0]; % green
textcolorfast = [255 200 0]; % red
textcolorslow = [255 0 0]; % red

%----------------------------------------------------------------------
%                       Response device settings (Keyboard)
%----------------------------------------------------------------------

% Response buttons (to produce sound)
Sound1 = KbName('left_shift');
Sound2 = KbName('left_control');
[keyIsDown, secs, keyCode] = KbCheck;
keyCode = zeros(1,size(keyCode,2));
SoundButtonTypes_LogicalArray = keyCode; SoundButtonTypes_LogicalArray(Sound1) = 1; SoundButtonTypes_LogicalArray(Sound2) = 1;
SoundButtonTypes_LogicalArray = logical(SoundButtonTypes_LogicalArray);

if counterbalancing == 1 % then uppercase button produces low sound
    SoundButtonTypes = [Sound1; Sound2];
elseif counterbalancing == 2 % then uppercase button produces high sound
    SoundButtonTypes = [Sound2; Sound1];
end

% Response buttons (for Judgement of Agency)
yes = KbName('4');
halfYes = KbName('7');
no = KbName('5');
halfNo= KbName('8');
keyCode = zeros(1,size(keyCode,2));
JoAButtonTypes_LogicalArray = keyCode;
JoAButtonTypes_LogicalArray(yes) = 1; JoAButtonTypes_LogicalArray(halfYes) = 1; JoAButtonTypes_LogicalArray(no) = 1; JoAButtonTypes_LogicalArray(halfNo) = 1;
JoAButtonTypes_LogicalArray = logical(JoAButtonTypes_LogicalArray);

JoAButtonTypes = [halfYes yes halfNo no];

% To start next trial
startTrial = KbName('right'); % flecha derecha


%----------------------------------------------------------------------
%                       Sound settings
%----------------------------------------------------------------------
FS=44100; % Sampling rate for sounds
headphone_factorH = 0.648;
headphone_factorL = 0.898;
starting_dB = 70;
intfactorH = 10^((starting_dB-100)/20)*headphone_factorH;  
intfactorL = 10^((starting_dB-100)/20)*headphone_factorL;  

% Load the sound and resample (needed for use in the lab)
SH = wavread ([Sounds_path '\SH_IA.wav']);
SL = wavread ([Sounds_path '\SL_IA.wav']);

SH = SH*intfactorH;
SL = SL*intfactorL;

% Start audio device
InitializePsychSound; % Initialize the sound device
paHandle = PsychPortAudio('Open',[],1,3,FS,2, [], [], [6 7]); % Open the Audio port and get a handle to refer to it in subsequent calls
% [6 7] are the headphones channels in the lab
Soundbuffers(1) = PsychPortAudio('Createbuffer',paHandle,SL);
Soundbuffers(2) = PsychPortAudio('Createbuffer',paHandle,SH);

%----------------------------------------------------------------------
%                       Prepare trials
%----------------------------------------------------------------------
nTrials = 20; % Trials per block


%% START EXPERIMENT

%----------------------------------------------------------------------
%                       Display initial instructions
%----------------------------------------------------------------------

[EXPWIN] = Screen('OpenWindow',0,backcolor); % Open fullscreen
Instr = ['Aprieta flecha derecha para empezar...' ];
DrawFormattedText(EXPWIN, Instr, 'center', 'center');
DrawFormattedText(EXPWIN, Instr, 'center', 'center');
Screen('Flip', EXPWIN);

% Start by pressing 0
keyCode = zeros(1,size(keyCode,2));
while not(keyCode(startTrial))
    [keyIsDown, secs, keyCode] = KbCheck;
end

% ----------------------------------------------------------------------
%                        Trial Loop preparation
% -----------------------------------------------------------------------

% Initialize loop variables
% Info we want to collect in Log to give feedback: 1) which buton was pressed, 2) what was the RT
Log = zeros (nTrials,2); % column 1: button (1 = left, 2 = right); column 2: cue-to-press RT
TrialEndTime = GetSecs; % dummy for first trial
SoundTime = []; % No sound played yet
Button = []; % No button pressed yet
JoA = []; % No JoA given yet
iTTrial = 1; % set test trial counter to start reading from TTrialArray
TrialTypeArray = 1; %Acquisition

% ----------------------------------------------------------------------
%                        TRIAL LOOP
% -----------------------------------------------------------------------

for iTrial = 1:nTrials
        
    % fixation has not changed color yet
    redflag = 0;
    greenflag = 0;
    orangeflag = 0;
    
    % flush kepresses
    keyCode = zeros(1,size(keyCode,2));
        
   % all ADQUISITION TRIALS
        % Show fixation cross (CUE)
        Screen('TextSize', EXPWIN, 55); % set text size
        Screen('TextFont', EXPWIN, 'Helvetica');
        DrawFormattedText(EXPWIN, '+', 'center', 'center',textcolor);
        FixTimeON = Screen('Flip',EXPWIN,TrialEndTime+(randsample(ITI,1))); % Draw fixation cross after waiting a random ITI from the end of previous trial

        %%%%%%%%%%%%%%%%%%%%%% PRESS WINDOW START %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % while we are in the response window, we continuously check for a button press to
        % present the sound and give feedback on the timing
        
        while GetSecs < (FixTimeON + RespWinDur)
            
            % if not press was detected yet, check for button press
            if isempty(Button)
                [keyIsDown, secs, keyCode] = KbCheck;
                
                % IF the button is pressed...
                if sum(keyCode(SoundButtonTypes_LogicalArray))
                    
                    % get the time and find out which button
                    PressTime = GetSecs;
                    Button = find(keyCode);
                    Button = intersect(Button,SoundButtonTypes);
                    % present the corresponding sound
                    Sound = find(Button == SoundButtonTypes); % 1 (low) or 2 (high), depending on counterbalancing
                    PsychPortAudio('FillBuffer', paHandle,Soundbuffers(Sound)); % Fill the audio buffer with sound corresponding to the button pressed
                    SoundTime = PsychPortAudio('Start', paHandle,1,1,1); % Play sound immediately
                    
                    % and show feedback (change fixation color)
                    if (PressTime - FixTimeON) < RTlimits(1) % too fast
                        feedbackcolor = textcolorfast;
                        orangeflag = 1;
                    elseif (PressTime - FixTimeON) > RTlimits(2) % too slow
                        feedbackcolor = textcolorslow;
                        redflag = 1;
                    else
                        feedbackcolor = textcolorcorrect; % on time
                        greenflag = 1;
                    end
                    DrawFormattedText(EXPWIN, '+', 'center', 'center',feedbackcolor);
                    Screen('Flip',EXPWIN,0);
                    
                end
            end
            
            % if RTlimits window has ended, and no press was recorded, paint the fixation red until end of
            % response window
            if (GetSecs > (FixTimeON + RTlimits(2))) && redflag == 0 && greenflag == 0 && orangeflag == 0
                DrawFormattedText(EXPWIN, '+', 'center', 'center',textcolorslow);
                Screen('Flip',EXPWIN,0);
                redflag = 1;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%% PRESS WINDOW END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Erase the fixation cross until next trial starts
        FixTimeOFF = Screen('Flip',EXPWIN,0); % Fix cross dissappears at FixTimeON + RespWinDur
        
        % if there was no press within the window, code for this as a miss:
        if isempty(Button)
            Button = 0;
            PressTime = NaN;
            Sound = 0;
            SoundTime = NaN;
        end
        
        % Set variables for log
        TrialEndTime = FixTimeON + Fix2QTime;
        JoA = NaN;
        JoATime = NaN;
        JoART = NaN;
        Condition = 1;
        Order = 1;
        RelSoundTime = SoundTime - PressTime;
        IntendedBin = NaN;
        ActualBin = NaN;
        QTime = NaN;
         
    % ----------------------------------------------------------------------
    %                         LOGFILE
    % -----------------------------------------------------------------------

    % Log press (or no press) in Log array for feedback at the end
    Log(iTrial,2)= PressTime - FixTimeON; % cue to press RT
    Log(iTrial,1)= Button; % which button was pressed
    
    % Write trial in logfile
    % Prepare string
    TrialLogString = [sprintf('%d',iTrial) ... % Column 1: Trial number
        sprintf('\t%02d', Subj) ... % Column 2: Subject number
        sprintf('\t%d', TrialTypeArray)... % Column 3: Trial Type (1 = acquisition, 2 = test)
        sprintf('\t%02d', Block) ... % Column 4: Block number (1-16)
        sprintf('\t%02d', Repetition) ... % Column 5: Repetition number
        sprintf('\t%5.7f', FixTimeON)... % Column 6: Cue time (resp window start)
        sprintf('\t%5.7f', FixTimeOFF)... % Column 7: Cue turns off (resp window end)
        sprintf('\t%d', Button) ... % Column 8: Button pressed
        sprintf('\t%5.7f', PressTime)... % Column 9: Button press time
        sprintf('\t%d', Sound)... % Column 10: Presented sound (1 = low, 2 = high)
        sprintf('\t%5.7f', SoundTime) ... % Column 11: Sound time
        sprintf('\t%1.3f', Log(iTrial,2))... % Column 12: cue to press RT
        sprintf('\t%5.7f', QTime) ... % Column 13: Question start time
        sprintf('\t%d', JoA)... % Column 14: Agency judgement (1-4)
        sprintf('\t%5.7f', JoATime) ... % Column 15: Agency judgement time
        sprintf('\t%1.3f', JoART)... % Column 16: Question to JoA RT
        sprintf('\t%d', Condition)... % Column 17: Condition (1 = congruent, 0 = incongruent, 2 = no press)
        sprintf('\t%d', Order)... % Column 18: Sound-Button order (-1 = sound first; 1 = button first)
        sprintf('\t%1.3f', RelSoundTime)... % Column 19: time of sound relative to time of press (-100 to +60)
        sprintf('\t%d', IntendedBin)... % Column 20: Time bin we aimed at (1-8)
        sprintf('\t%d', ActualBin)];% Column 21: time bin we hit (1-8)
    
    % Write string to logfile as a new line
    fprintf(LOGFILE,'\n%s',TrialLogString);
    
    % Clear loop vars
    Button = [];
    PressTime = [];
    Sound = [];
    JoA = [];
    SoundTime = [];
end

%% BLOCK FEEDBACK

%Define variables for feedback
LogButton = Log(:,1);
LogRT = Log(:,2);

% calculate and display percentage uppercase and control buttons
Button1percent = round((size(LogButton(LogButton==160),1)/size(LogButton,1))*100);
Button2percent = round((size(LogButton(LogButton==162),1)/size(LogButton,1))*100);

% Calculate and display mean RT and times too fast and times too slow
meanRT = nanmean(LogRT);
toofast = size(LogRT(LogRT<RTlimits(1)),1);
tooslow = size(LogRT(LogRT>RTlimits(2)),1);

% Reduce text size for text
Screen('TextSize', EXPWIN, 20); % set text size
Screen('TextFont', EXPWIN, 'Helvetica');
% Display feedback to the subject
Instr3 = ['media TR = ' sprintf('%1.1f',meanRT) '\n\n\ndemasiado lento = ' num2str(tooslow) '      demasiado rapido = ' num2str(toofast) '\n\n\nBoton1 = ' num2str(Button1percent) '%   Boton2 = ' num2str(Button2percent) '%'];
DrawFormattedText(EXPWIN, Instr3, 'center', 'center',textcolor);
DrawFormattedText(EXPWIN, Instr3, 'center', 'center',textcolor);
Screen('Flip',EXPWIN,0) % StartTime = 0 for now

WaitSecs(5);
disp('Press any key to finish');
pause;

%% Close things
fclose(LOGFILE);
PsychPortAudio('Close'); % Close the audio device
Screen('Close',EXPWIN);

end

