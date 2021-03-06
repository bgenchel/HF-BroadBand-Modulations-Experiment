% Runs subject through the NEO 5 point personality test
% Bought Oct, 2005
% 

allscore = zeros(30,8); % initialize matrix of correct size
p1 = 1; % matrix place indicator for row (increment first 1:30, then start over at pl2...)
p2 = 1; % matrix place indicator for col

quests = {'I am not a worrier.'; ...
         'I really like most people I meet.'; ...
        'I have a very active imagination.';...
        'I tend to be cynical and skeptical of others'' intentions.';...
        'I''m known for my prudence and common sense.';...
        'I often get angry at the way people treat me.';...
        'I shy away from crowds of people.';...
        'Aesthetic and artistic concerns aren''t very important to me.';...
        'I''m not crafty or sly.';...
        'I would rather keep my options open than plan everything in advance.';...
        'I rarely feel lonely or blue.';...
        'I am dominant, forceful, and assertive.';...
        'Without strong emotions, life wold be uninteresting to me.';...
        'Some people think I''m selfish and egotistical.';...
        'I try to perform all the tasks assigned to me conscientiously.';...
        'In dealing with other people, I always dread making a social blunder.';...
        'I have a leisurely style in work and play.';...
        'I''m pretty set in my ways.';...
        'I would rather cooperate with others than compete with them.';...
        'I am easy-going and lackadaisical.';...
        'I rarely overindulge in anything.';...
        'I often crave excitement.';...
        'I often enjoy playing with theories or abstract ideas.';...
        'I don''t mind bragging about my talents and accomplishments.';...
        'I''m pretty good about pacing myself so as to get things done on time.';...
        'I often feel helpless and want someone else to solve my problems.';...
        'I have never literally jumped for joy.';...
        'I believe letting students hear controversial speakers can only confuse and mislead them.';...
        'Political leaders need to be more aware of the human side of their policies.';...
        'Over the years I''ve done some pretty stupid things.';...
        'I am easily frightened.';...
        'I don''t get much pleasure from chatting with people.';... 
        'I try to keep all my thoughts directed along realistic lines and avoid flights of fancy.';...
        'I believe that most people are basically well-intentioned.';...
        'I don''t take civic duties like voting very seriously.';...
        'I''m an even-tempered person.';...
        'I like to have a lot of people around me.';...
        'I am sometimes completely abosorbed in music I am listening to.';...
        'If necessary, I am willing to manipulate people to get what I want.';...
        'I keep my belongings neat and clean.';...
        'Sometimes I feel completely worthless.';...
        'I sometimes fail to assert myself as much as I should.';...
        'I rarely experience strong emotions.';...
        'I try to be courteous to everyone I meet.';...
        'Sometimes I''m not as dependable or reliable as I should be.';...
        'I seldom feel self-conscious when I''m around people.';...
        'When I do things, I do them vigorously.';...
        'I think it''s interesting to learn and develop new hobbies.';...
        'I can be sarcastic and cutting when I need to be.';...
        'I have a clear set of goals and work toward them in an orderly fashion.';...
        'I have trouble resisting my cravings.';...
        'I wouldn''t enjoy vacationing in Las Vegas.';...
        'I find philosophical arguments boring.';...
        'I''d rather not talk about myself and my achievements.';...
        'I waste a lot of time before settling down to work.';...
        'I feel I am capable of coping with most of my problems.';...
        'I have sometimes experienced intense joy or ecstasy.';...
        'I believe that laws and social policies should change to reflect the needs of a changing world.';...
        'I''m hard-headed and tough-minded in my attitudes.';...
        'I think things through before coming to a decision.';...
        'I rarely feel fearful or anxious.';...
        'I''m known as a warm and friendly person.';...
        'I have an active fantasy life.';...
        'I believe that most people will take advantage of you if you let them.';...
        'I keep myself informed and ususally make intelligent decisions.';...
        'I am known as hot-blooded and quick-tempered.';...
        'I usually prefer to do things alone.';...
        'Watching ballet or modern dance bores me.';...
        'I couldn''t deceive anyone even if I wanted to.';...
        'I am not a very methodical person.';...
        'I am seldom sad or depressed.';...
        'I have often been a leader of groups I have belonged to.';...
        'How I feel about things is important to me.';...
        'Some people think of me as cold and calculating.';...
        'I pay my debts promptly and in full.';...
        'At times I have been so ashamed I just want to hide.';...
        'My work is likely to be slow but steady.';...
        'Once I find the right way to do something, I stick to it.';...
        'I hesitate to express my anger even when it''s justified.';...
        'When I start a self-improvement program, I usually let it slide after a few days.';...
        'I have little difficulty resisting temptation.';...
        'I have sometimes done things just for "kicks" or "thrills".';...
        'I enjoy solving problems or puzzles.';...
        'I''m better than most people, and I know it.';...
        'I am a productive person who always gets the job done.';...
        'When I''m under a great deal of stress, sometimes I feel like I''m going to pieces.';...
        'I am not a cheerful optimist.';...
        'I believe we should look to our religious authorities for decisions on moral issues.';...
        'We can never do too much for the poor and elderly.';...
        'Occasionally I act first and think later.';...
        'I often feel tense and jittery.';...
        'Many people think of me as somewhat cold and distant.';...
        'I don''t like to waste my time daydreaming.';...
        'I think most of the people I deal with are honest and trustworthy.';...
        'I often come into situations without being fully prepared.';...
        'I am not considered a touchy or temperamental person.';...
        'I really feel the need for other people if I am by myself for long.';...
        'I am intrigued by the patterns I find in art and nature.';...
        'Being perfectly honest is a bad way to do business.';...
        'I like to keep everything in its place so I know just where it is.';...
        'I have sometimes experienced a deep sense of guilt or sinfulness.';...
        'In meetings, I usually let others do the talking.';...
        'I seldom pay much attention to my feelings of the moment.';...
        'I generally try to be thoughtful and considerate.';...
        'Sometimes I cheat when I play solitaire.';...
        'It doesn''t embarrass me too much if people ridicule and tease me.';...
        'I often feel as if I''m bursting with energy.';...
        'I often try new and foreign foods.';...
        'If I don''t like people, I let them know it.';...
        'I work hard to accomplish my goals.';...
        'When I am having my favorite foods, I tend to eat too much.';...
        'I tend to avoid movies that are shocking or scary.';...
        'I sometimes lose interest when people talk about very abstract, theoretical matters.';...
        'I try to be humble.';...
        'I have trouble making myself do what I should.';...
        'I keep a cool head in emergencies.';...
        'Sometimes I bubble with happiness.';...
        'I believe that the different ideas of right and wrong that people in other societies have may be valid for them.';...
        'I have no sympathy for panhandlers.';...
        'I always consider the consequences before I take action.';...
        'I''m seldom apprehensive about the future.';...
        'I really enjoy talking to people.';...
        'I enjoy concentrating on a fantasy or daydream and exploring all its possibilities, letting it grow and develop.';...
        'I''m suspicious when someone does something nice for me.';...
        'I pride myself on my sound judgment.';...
        'I often get disgusted with people I have to deal with.';...
        'I prefer jobs that let me work alone without being bothered be other people.';...
        'Poetry has little or no effect on me.';...
        'I would hate to be thought of as a hypocrite.';...
        'I never seem able to get organized.';...
        'I tend to blame myself when anything goes wrong.';...
        'Other people often look to me to make decisions.';...
        'I experience a wide range of emotions or feelings.';...
        'I''m not known for my generosity.';...
        'When I make a commitment, I can always be counted on to follow through.';...
        'I often feel inferior to others.';...
        'I''m not as quick and lively as other people.';...
        'I prefer to spend my time in familiar surroundings.';...
        'When I''ve been insulted, I just try to forgive and forget.';...
        'I don''t feel like I''m driven to get ahead.';...
        'I seldom give in to my impulses.';...
        'I like to be where the action is.';...
        'I enjoy working on "mind-twister"-type puzzles.';...
        'I have a very high opinion of myself.';...
        'Once I start a project, I almost always finish it.';...
        'It''s often hard for me to make up my mind.';...
        'I don''t consider myself especially "light-hearted."';...
        'I believe that loyalty to one''s ideals and principles is more important than "open-mindedness."';...
        'Human need should always take priority over economic situations.';...
        'I often do things on the spur of the moment.';...
        'I often worry about things that might go wrong.';...
        'I find it easy to smile and be outgoing with strangers.';...
        'If I feel my mind starting to drift off into daydreams, I usually get busy and start concentrating on some work or activity instead.';...
        'My first reaction is to trust people.';...
        'I don''t seem to be completely successful at anything.';...
        'It takes a lot to get me mad.';...
        'I''d rather vacation at a popular beach than as isolated cabin in the woods.';...
        'Certain kinds of music have an endless fascination for me.';...
        'Sometimes I trick people into doing what I want.';...
        'I tend to be somewhat fastidious or exacting.';...
        'I have a low opiniong of myself.';...
        'I would rather go my own way than be a leader of others.';...
        'I seldom notice the moods or feelings that different environments produce.';...
        'Most people I know like me.';...
        'I adhere strictly to my ethical principles.';...
        'I feel comfortable in the presence of my bosses or other authorities.';...
        'I usually seem to be in a hurry.';...
        'Sometimes I make changes around the house just to try something different.';...
        'If someone starts a fight, I''m ready to fight back.';...
        'I strive to achieve all I can.';...
        'I sometimes eat myself sick.';...
        'I love the excitement of roller coasters.';...
        'I have little interest in speculating on the nature of the universe or the human condition.';...
        'I feel that I am no better than others, no matter what their condition.';...
        'When a project gets too difficult, I''m inclined to start a new one.';...
        'I can handle myself pretty well in a crisis.';...
        'I am a cheerful, high-spirited person.';...
        'I consider myself broad-minded and tolerant of other people''s lifestyles.';...
        'I believe all human beings are worthy of respect.';...
        'I rarely make hasty decisions.';...
        'I have fewer fears than most people.';...
        'I have strong emotional attachments to my friends.';...
        'As a child I rarely enjoyed games of make believe.';...
        'I tend to assume the best about people.';...
        'I''m a very competent person.';...
        'At times I have felt bitter and resentful.';...
        'Social gatherings are usually boring to me.';...
        'Sometimes when I am reading poetry or looking at a work of art, I feel a chill or wave of excitement.';...
        'At time I bully or flatter people into doing what I want them to.';...
        'I''m not compulsive about cleaning.';...
        'Sometimes things look pretty bleak and hopeless to me.';...
        'In conversations, I tend to do most of the talking.';...
        'I find it easy to empathize - to feel myself what others are feeling.';...
        'I think of myself as a charitable person.';...
        'I try to do jobs carefully, so they won''t have to be done again.';...
        'If I have said or done the wrong thing to someone, I can hardly bear to face them again.';...
        'My life is fast-paced.';...
        'On a vacation, I prefer going back to a tried and true spot.';...
        'I''m hard-headed and stubborn.';...
        'I strive for excellence in everything I do.';...
        'Sometimes I do things on impulse that I later regret.';...
        'I''m attracted to bright colors and flashy styles.';...
        'I have a lot of intellectual curiosity.';...
        'I would rather praise others than be praised myself.';...
        'There are so many little jobs that need to be done that I sometimes just ignore them all.';...
        'When everything seems to be going wrong, I can still make good decisions.';...
        'I rarely use words like "fantastic!" or "sensational!" to describe my experiences.';...
        'I think that if people don''t know what they believe in by the time they''re 25, there''ssomething wrong with them.';...
        'I have sympathy for others less fortunate than me.';...
        'I plan ahead carefully when I go on a trip.';...
        'Frightening thoughts sometimes come into my head.';...
        'I take a personal interest in the people I work with.';...
        'I would have difficulty just letting my mind wander without control or guidance.';...
        'I have a good deal of faith in human nature.';...
        'I am efficient and effective at my work.';...
        'Even minor annoyances can be frustrating to me.';...
        'I enjoy parties with lots of people.';...
        'I enjoy reading poetry that emphasizes feelings and images more than story lines.';...
        'I pride myself on my shrewdness in handling people.';...
        'I spend a lot of time looking for things I''ve misplaced.';...
        'Too often, when things go wrong, I get discouraged and feel like giving up.';...
        'I don''t find it easy to take charge of a situation.';...
        'Odd things - like certain scents or the names of distant places - can evoke strong moods in me.';...
        'I go out of my way to help others if I can.';...
        'I''d really have to be sick before I''d miss a day of work.';...
        'When people I know do foolish things, I get embarrassed for them.';...
        'I am a very active person.';...
        'I follow the same route when I go someplace.';...
        'I often get into arguments with my family and co-workers.';...
        'I''m something of a "workaholic."';...
        'I am always able to keep my feelings under control.';...
        'I like being part of the crowd at sporting events.';...
        'I have a wide range of intellectual interests.';...
        'I''m a superior person.';...
        'I have a lot of self-discipline.';...
        'I''m pretty stable emotionally.';...
        'I laugh easily.';...
        'I believe that the "new morality" of permissiveness in no morality at all.';...
        'I would rather be known as "merciful" than as "just."';...
        'I think twice before I answer a question.'};
quests = {'I am not a worrier.'; ...
         'I believe that laws and social policies should change to reflect the needs of a changing world.'};
          
clear qidx          
qidx = ones(1,length(quests));% initialize to all 4:0 questions
qidx([2,3,5,6,9,12,13,15,16,19,22,23,25,26,29,31,34,37,38,40,41,44,47,48,50,51,54,57,58,60,62,63,65,66,69,72,73,75,76,79,82,83,85,86,89,91,94,97,98,100,101,104,107,108,110,111,114,117,118,120,122,123,125,126,129,131,132,133,135,136,139,142,143,145,146,149,151,152,154,157,158,160,161,164,165,167,168,170,171,172,174,177,178,179,180,182,184,185,186,188,191,192,193,194,195,196,197,200,201,202,203,204,209,210,211,212,214,215,216,217,218,221,223,224,225,226,227,230,232,233,235,237,239,240]) = 2;    % make 0:4 quests 2

sdvals = [4,0];
dvals =  [3,1];
avals =  [1,3];
savals = [0,4];

h_fig = figure('MenuBar','none');
set(gcf,'color',[.5 .8 1]); set(gcf,'position',[320 390 860 420]);set(gcf,'units','normalized');


for q = 1:length(quests)
    callfunc1 = ['allscore(',int2str(p1),',',int2str(p2),') = get(h_box1,''userdata'');set(txt_box,''visible'',''off'');'];
    callfunc2 = ['allscore(',int2str(p1),',',int2str(p2),') = get(h_box2,''userdata'');set(txt_box,''visible'',''off'');'];
    callfunc3 = ['allscore(',int2str(p1),',',int2str(p2),') = get(h_box3,''userdata'');set(txt_box,''visible'',''off'');'];
    callfunc4 = ['allscore(',int2str(p1),',',int2str(p2),') = get(h_box4,''userdata'');set(txt_box,''visible'',''off'');'];
    callfunc5 = ['allscore(',int2str(p1),',',int2str(p2),') = get(h_box5,''userdata'');set(txt_box,''visible'',''off'');'];
    txt_box = uicontrol(h_fig,'units','normalized','Position',[0 0 .6 .3],'style','text','FontSize',18,'string',quests{q},'visible','off');
    newpos = get(txt_box,'extent') ;
    newpos = newpos(3:4)*1.1;
    txt_box = uicontrol(h_fig,'units','normalized','Position',[.5-(newpos(1)/2) .6 newpos(1) newpos(2)],'style','text','FontSize',18,'string',quests{q},'backgroundcolor',[.5 .8 1]);

    sdval = sdvals(qidx(q));
    dval =  dvals(qidx(q));
    aval =  avals(qidx(q));
    saval = savals(qidx(q));

    h_box1 = uicontrol(h_fig,'units','normalized','Position',[.02 .2 .17 .07],'string','Strongly Disagree','style','pushbutton','userdata',sdval,'ForegroundColor','b','FontSize',12,'backgroundcolor',[.6 .6 .6],'callback',callfunc1);

    h_box2 = uicontrol(h_fig,'units','normalized','Position',[.22 .2 .17 .07],'string','Disagree','style','pushbutton','userdata',dval,'ForegroundColor','b','FontSize',12,'backgroundcolor',[.6 .6 .6],'callback',callfunc2);

    h_box3 = uicontrol(h_fig,'units','normalized','Position',[.42 .2 .17 .07],'string','Neutral','style','pushbutton','userdata',[2],'ForegroundColor','b','FontSize',12,'backgroundcolor',[.6 .6 .6],'callback',callfunc3);

    h_box4 = uicontrol(h_fig,'units','normalized','Position',[.62 .2 .17 .07],'string','Agree','style','pushbutton','userdata',aval,'ForegroundColor','b','FontSize',12,'backgroundcolor',[.6 .6 .6],'callback',callfunc4);

    h_box5 = uicontrol(h_fig,'units','normalized','Position',[.82 .2 .17 .07],'string','Strongly Agree','style','pushbutton','userdata',saval,'ForegroundColor','b','FontSize',12,'backgroundcolor',[.6 .6 .6],'callback',callfunc5);

    waitfor(txt_box,'visible','off');
    
    p1 = p1+1;
    if p1 == 31|p1==61|p1==91|p1==121|p1==151|p1==181|p1==211
        p2 = p2+1;
        p1 = 1;
    end;   
    clf
end;

close

categ = {'Anxiety','Angry Hostility','Depression','Self-Consciousness','Impulsiveness','Vulnerability','Warmth','Gregariousness','Assertiveness','Activity','Excitement-Seeking','Positive Emotions','Fantasy','Aesthetics','Feelings','Actions','Ideas','Values', 'Trust','Straightforwardness','Altruism','Compliance','Modesty','Tender-Mindedness', 'Competence','Order','Dutifulness  ','Achievement Striving','Self-Discipline','Deliberation'};

ranges = zeros(30,4);
ranges(1:30,1) = [24,21,22,22,24,18,30,25,23,25,24,28,24,27,28,23,26,27,28,29,30,26,26,26,28,26,30,26,29,24]';
ranges(1:30,2) = [19,16,16,18,19,14,26,20,19,21,19,24,19,22,24,19,21,23,24,25,27,22,22,23,24,22,2,22,25,20]';
ranges(1:30,3) = [13,10,10,13,14,9,22,15,13,16,13,19,14,16,19,15,16,19,20,20,23,18,18,20,20,17,22,18,20,15]';
ranges(1:30,4) = [8,6,5,9,10,5,18,10,8,11,8,14,9,11,15,12,11,15,16,16,20,14,14,17,17,13,18,14,15,11]';

for rw = 1:size(ranges,1)
    if rw < 7
        N(1,rw) = sum(allscore(rw,:));
    end;
    if rw > 6 & rw < 13
        E(1,rw-6) = sum(allscore(rw,:));
    end;
    if rw > 12 & rw < 19
        O(1,rw-12) = sum(allscore(rw,:));
    end;
    if rw > 18 & rw < 25
        A(1,rw-18) = sum(allscore(rw,:));
    end;
    if rw > 24
        C(1,rw-24) = sum(allscore(rw,:));
    end;

    if sum(allscore(rw,:)) >= ranges(rw,1)    
        msg = 'very high';
    elseif sum(allscore(rw,:)) >= ranges(rw,2) 
        msg = 'high';
    elseif sum(allscore(rw,:)) >= ranges(rw,3) 
        msg = 'average';
    elseif sum(allscore(rw,:)) >= ranges(rw,4) 
        msg = 'low';
    else
        msg = 'very low';
    end;
if length(categ{rw}) > 16
        fprintf('\n%s:\t%s',categ{rw},msg);
elseif length(categ{rw}) < 8
    fprintf('\n%s\t\t:\t%s',categ{rw},msg);
else
    fprintf('\n%s\t:\t%s',categ{rw},msg);    
end;
end;
fprintf('\n');
totranges(1:5,1) = [117,139,138,151,151]';
totranges(1:5,2) = [96,121,121,137,133]';
totranges(1:5,3) = [72,101,102,121,113]';
totranges(1:5,4) = [50,82,85,107,96]';

N(1,end+1) = sum(N(1:6));
E(1,end+1) = sum(E(1:6));
O(1,end+1) = sum(O(1:6));
A(1,end+1) = sum(A(1:6));
C(1,end+1) = sum(C(1:6));

totcats = {'Neuroticism','Extroversion','Openness  ','Agreeableness','Conscientiousness'};
alltots = [N(end),E(end),O(end),A(end),C(end)];

for tt = 1:length(alltots)
    if alltots(tt) >= totranges(tt,1)    
        msg = 'very high';
    elseif alltots(tt) >= ranges(tt,2) 
        msg = 'high';
    elseif alltots(tt) >= ranges(tt,3) 
        msg = 'average';
    elseif alltots(tt) >= ranges(tt,4) 
        msg = 'low';
    else
        msg = 'very low';
    end;
    if length(totcats{tt}) > 16
        fprintf('\n%s:\t%s',totcats{tt},msg);
    else        
        fprintf('\n%s\t:\t%s',totcats{tt},msg);
    end;
end;
