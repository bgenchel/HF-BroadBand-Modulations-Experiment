% runs through datasets and collects emotion intensity self-ratings
% eeglab
filnm = {'awe.set' 'awe2.set' 'frust.set' 'frust2.set' 'joy.set' 'joy2.set' 'anger.set' 'anger2.set' 'sad.set' 'sad2.set' 'surprise.set' 'surprise2.set' 'happy.set' 'happy2.set' 'fear.set' 'fear2.set' 'love.set' 'love2.set' 'jealousy.set' 'jealousy2.set' 'compassion.set' 'compassion2.set' 'content.set' 'grief.set' 'grief2.set' 'relief.set' 'relief2.set'};%'content2.set'
 clear cond rate
for p = 1:length(filnm)
    EEG = pop_loadset( filnm{p},'/data/common1/emotion/rr80/');
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    cond{1,p} = EEG.event(1).type;
    rate(p,1) = EEG.event(4).type;
    ALLEEG=[];
end;
for p = 1:length(rate)
    fprintf('Condition: %s\t', cond{1,p});
    fprintf('Self-Rating: %d\n',rate(p,1));
end;

ap80:

Condition: awe          Self-Rating: 6
Condition: awe 2        Self-Rating: 4
Condition: frust        Self-Rating: 7
Condition: frust 2      Self-Rating: 6
Condition: joy          Self-Rating: 6
Condition: joy 2        Self-Rating: 4
Condition: anger        Self-Rating: 7
Condition: anger 2      Self-Rating: 5
Condition: sad          Self-Rating: 5
Condition: sad 2        Self-Rating: 4
Condition: surprise     Self-Rating: 3
Condition: surprise 2   Self-Rating: 6
Condition: happy        Self-Rating: 3
Condition: happy 2      Self-Rating: 6
Condition: fear         Self-Rating: 5
Condition: fear 2       Self-Rating: 5
Condition: love         Self-Rating: 6
Condition: love 2       Self-Rating: 6
Condition: jealousy     Self-Rating: 6
Condition: jealousy 2   Self-Rating: 5
Condition: compassion   Self-Rating: 5
Condition: compassion 2 Self-Rating: 4
Condition: content      Self-Rating: 6
Condition: content 2    Self-Rating: 6
Condition: grief        Self-Rating: 4
Condition: grief 2      Self-Rating: 4
Condition: relief       Self-Rating: 6
Condition: relief 2     Self-Rating: 6

rr80:

Condition: awe          Self-Rating: 3
Condition: awe          Self-Rating: 6
Condition: frust        Self-Rating: 7
Condition: frust        Self-Rating: 7
Condition: joy          Self-Rating: 6
Condition: joy          Self-Rating: 6
Condition: anger        Self-Rating: 8
Condition: anger        Self-Rating: 8
Condition: sad          Self-Rating: 9
Condition: sad          Self-Rating: 7
Condition: surprise     Self-Rating: 1
Condition: surprise     Self-Rating: 1
Condition: happy        Self-Rating: 8
Condition: happy        Self-Rating: 8
Condition: fear         Self-Rating: 7
Condition: fear         Self-Rating: 4
Condition: love         Self-Rating: 8
Condition: love         Self-Rating: 7
Condition: jealousy     Self-Rating: 6
Condition: jealousy     Self-Rating: 4
Condition: compassion   Self-Rating: 5
Condition: compassion   Self-Rating: 2
Condition: content      Self-Rating: 8
Condition: grief        Self-Rating: 9
Condition: grief        Self-Rating: 1
Condition: relief       Self-Rating: 8
Condition: relief       Self-Rating: 8
