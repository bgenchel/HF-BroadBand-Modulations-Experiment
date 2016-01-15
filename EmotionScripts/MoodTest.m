% Runs a current mood test by having subjects click on successive figures with 5 point scales

savedat = '/data/common4/nf/';% if not [], saves data to this directory in subject subdir

h_fig = figure('MenuBar','none');
set(gcf,'color',[.5 .8 1]); set(gcf,'position',[320 390 860 520]);set(gcf,'units','normalized');

ylev = .08;  lft = .1;
txt_box = uicontrol(h_fig,'units','normalized','Position',[lft .5 .8 .4],'style','text','FontSize',18,'string','You will see a series of windows each presenting a single word that describes a particular feeling or emotion. Read each word and then click on the red line to indicate to what extent you feel this way RIGHT NOW, that is, at the present moment. Use the scale indicated in each figure (1: slightly or not at all, to 5: extremely) to rate your answers.','backgroundcolor',[.5 .8 1],'horizontalalignment','left');

ylev = .41;  
txt_box = uicontrol(h_fig,'units','normalized','Position',[lft ylev .4 .08],'style','text','FontSize',18,'string','Please type in your subject id: ','backgroundcolor',[.5 .8 1]);
id_box = uicontrol(h_fig,'units','normalized','Position',[lft+.45 ylev+.02 .15 .07],'style','edit','FontSize',20);

txt_box = uicontrol(h_fig,'units','normalized','Position',[lft ylev-.08 .7 .08],'style','text','FontSize',14,'string','(initials and year of birth: ie, John Smith born in 1981 would be ''js81'')','backgroundcolor',[.5 .8 1],'horizontalalignment','left');

txt_box = uicontrol(h_fig,'units','normalized','Position',[lft+.082 ylev-.17 .4 .08],'style','text','FontSize',18,'string','Gender (f or m): ','backgroundcolor',[.5 .8 1]);
gend_box = uicontrol(h_fig,'units','normalized','Position',[lft+.45 ylev-.17+.02 .15 .07],'style','edit','FontSize',22);

end_box = uicontrol(h_fig,'units','normalized','Position',[.5-.1 ylev-.35 .2 .15],'string','BEGIN','ForegroundColor','r','FontWeight','bold','FontSize',16,'style','toggle');

waitfor(end_box,'value',1);
subjid = get(id_box,'string');
subjgend = get(gend_box,'string');
close


fprintf('\nSubject: %s',subjid);

words20 = {'afraid','scared','nervous','jittery','irritable','hostile','guilty','ashamed','upset','distressed','active','alert','attentive','determined','enthusiastic','excited','inspired','interested','proud','strong'};
         
wordsall = {'afraid','scared','frightened','nervous','jittery','shaky','angry','hostile','irritable','scornful','disgusted','loathing','guilty','ashamed','blameworthy','angry at self','disgusted with self','dissatisfied with self','sad','blue','downhearted','alone','lonely','upset','distressed', 'happy','joyful','delighted','cheerful','excited','enthusiastic','lively','energetic','proud','strong','confident','bold','daring','fearless','alert','attentive','concentrating','determined','active','inspired','interested','shy','bashful','sheepish','timid','sleepy','tired','sluggish','drowsy','calm','relaxed','at ease','amazed','surprised','astonished'};


%wordskim = {'interested','distressed','excited','upset','strong','guilty','scared','hostile','enthusiastic','proud','competent','stupid','resourceful','efficient','shame','smart','worthless','boastful','fearful','frightened','impatient','worrying','irritable','alert','ashamed','inspired','nervous','determined','attentive','jittery','active','afraid','egotistic','effective','conceited','inadequate','confident','incompetent','self-centered','panicky','shaky','tense','timid'};

%words = wordsall; tst = 1; % choose which list to use.
words = words20; tst = 2; % choose which list to use.
%words = {'interested','distressed'};

randwords = [1:length(words)];
randwords = shuffle(randwords);

figure; clear wscore
for wd = 1:length(words)
    word = randwords(wd);

    ph=plot([1 5],[0 0],'r-','linewidth',10); hold on;
    set(ph,'linewidth',10);
    set(gca,'xlim',[.99 5.01]);set(gca,'ylim',[-.5 .5]);
    axis('off');set(gcf,'color','w');
    sub = round(length(words{word})/2)/10;

    ph = text(3-sub,.35,words{word});set(ph,'fontsize',18);
    set(ph,'color','b');

    plot([1 1],[-.1 .1],'r-','linewidth',2);
    plot([2 2],[-.1 .1],'r-','linewidth',2);
    plot([3 3],[-.1 .1],'r-','linewidth',2);
    plot([4 4],[-.1 .1],'r-','linewidth',2);
    plot([5 5],[-.1 .1],'r-','linewidth',2);
    fs = 16;    
    ph=text(.5,-.24,['very slightly'; 'or not at all']);set(ph,'fontsize',fs);
    ph=text(1.8,-.2,'a little');set(ph,'fontsize',fs);
    ph=text(2.65,-.24,'moderately');set(ph,'fontsize',fs);
    ph=text(3.65,-.2,'quite a bit');set(ph,'fontsize',fs);
    ph=text(4.7,-.24,'extremely');set(ph,'fontsize',fs);
    fs = 18;
    ph=text(.97,.16,'1');set(ph,'fontsize',fs);
    ph=text(1.97,.16,'2');set(ph,'fontsize',fs);
    ph=text(2.97,.16,'3');set(ph,'fontsize',fs);
    ph=text(3.97,.16,'4');set(ph,'fontsize',fs);
    ph=text(4.97,.16,'5');set(ph,'fontsize',fs);

    [wscore(word) y] = ginput(1);
    clf
end;
close
for wd =1:length(words)
fprintf('\n%s   \t%1.3g',words{wd},wscore(wd));
end;
fprintf('\n\n');
if tst == 1
    NA = sum(wscore(1:25)); %basic negative emos
    PA = sum(wscore(26:46)); %basic positive emos
    shy = sum(wscore(47:50));
    fatigue = sum(wscore(51:54));
    serenity = sum(wscore(55:57));
    suprise = sum(wscore(58:60));
elseif tst == 2
    NA = sum(wscore(1:10)); %basic negative emos   (~15.5 +- 5.5) m +- sd
    PA = sum(wscore(11:20)); %basic positive emos  (~27.15 +- 7.5) m +- sd
end;

%str = [savedat,subjpath,'/MoodScorePre.mat wscore words subjgend']; eval(str);


words2 = {'interested','young','terrible','discouraged','sad','good','alive','strong','active','miserable','free','healthy','unhappy','peaceful','whole','safe','gloomy','merry','suffering','tormented','destroyed','clean','forlorn','wilted','hopeless','lucky','rejected','low','gay','inspired','awful','glad','lonely','lost','blue','enthusiastic','sunk','alone','fine','fit'};

words2={'interested' , 'terrible' ,   'discouraged'  ,  'sad',    'good', 'alive',    'strong',    'active',    'miserable',    'free',    'healthy',  'unhappy',    'peaceful',    'whole',    'safe',    'gloomy',    'merry',  'destroyed',   'hopeless',    'lucky',    'rejected',    'low',    'inspired',    'awful' ,   'glad',    'lost',    'blue',   'enthusiastic',    'alone',    'fine',    'fit'};

words2 = {'afraid','nervous','tense','irritable','dissatisfied with self','sad','downhearted','lonely','distressed', 'happy','cheerful','enthusiastic','lively','energetic','proud','alert','inspired','interested','timid','tired','drowsy','calm','relaxed','at ease','amazed'};

h_fig = figure('MenuBar','none');
set(gcf,'color',[.5 .8 1]); set(gcf,'position',[320 390 960 520]);set(gcf,'units','normalized');
%statchk = ['clkword = get(h_box,''value'');'];
ph=text(.11,1.02,'Please click on all the adjectives that indicate how you feel right now:');axis('off');
set(ph,'fontsize',16);
randwords = [1:length(words2)];
randwords = shuffle(randwords);
ht = .06; incfac = .07;y = .88;
for wd = 1:length(words2)
    word = randwords(wd);
    if wd < ceil(length(words2)/3) + 1
        x = .11;
    elseif wd >= ceil(length(words2)/3)+1  & wd < ceil(length(words2)/3)*2+1
        x=.41;
    else
        x = .71;
    end;
    if wd == ceil(length(words2)/3)+1 | wd == ceil(length(words2)/3)*2 +1
        endbox = y; y = .88;    
    end;
    y = y - incfac;    
    h_box = uicontrol(h_fig,'units','normalized','Position',[x y .19 ht],'string',words2{word},'style','checkbox','value',[0],'ForegroundColor','b','Tag',words2{word},'FontSize',12);
    allbox{word} = h_box;
end;
end_box = uicontrol(h_fig,'units','normalized','Position',[x endbox-.05 .19 .07],'string','Finished','ForegroundColor','r','FontWeight','bold','FontSize',12,'style','toggle');
waitfor(end_box,'value',1);
for word = 1:length(words2)
    wordvals(word) = get(allbox{word},'value');
    if wordvals(word) == 1
        words2{word};
    end;    
end;
close  
fprintf('\nWords that describe how you feel right now:\n');
for wd =1:length(words2)
    if wordvals(wd) == 1
        fprintf('\n%s',words2{wd});
    end;
end;
fprintf('\n\n');

%str = ['save ',subjpath,'MoodScorePre.mat wscore words wordvals words2']; eval(str);
 
