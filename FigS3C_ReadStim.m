% ReadSpectraVS
close all
clear all

path='E:\'; 
pathsig=[path,'OptoInh\OutputSignals\'];
pathstim=[path,'OptoInh\STIMs\']; mkdir(pathstim)


%%%% SD
% event='est1'
% maxep=10800;
% mousenames=char('Ar1');days=['120421'];
% mousenames=char('Ar2');days=['080421'];
% mousenames=char('Ar4');days=['120421'];
% mousenames=char('Ar5');days=['060521'];
% mousenames=char('Ar8');days=['060521'];
% mousenames=char('Ar9');days=['100521'];
% mousenames=char('Gf1');days=['080421'];
% mousenames=char('Gf2');days=['080421'];
% mousenames=char('Gf4');days=['100521'];
% mousenames=char('Gf5');days=['060521'];

%%%% 1h continous stim
% event='est2'
% maxep=10800;
% mousenames=char('Gf1');days=['170421'];
event='est1';
maxep=10800;
mousenames=char('Ar1');days=['170421'];
% event='est1'
% maxep=10800;
% mousenames=char('Ar1');days=['170421D4'];

% %%%% 24hr stim
% event='est4'
% maxep=21600;
% mousenames=char('Ar1');days=['160421'];
% mousenames=char('Ar2');days=['150421'];
% mousenames=char('Ar4');days=['160421'];
% mousenames=char('Ar5');days=['140521'];
% mousenames=char('Ar8');days=['140521'];
% mousenames=char('Ar9');days=['150521'];
% mousenames=char('Gf1');days=['150421'];
% mousenames=char('Gf2');days=['150421'];
% mousenames=char('Gf4');days=['150521'];
% mousenames=char('Gf5');days=['140521'];


%%%%%% Wake Enhancement %50 times stim between ZT2 and ZT4
%%%% LPO
% mousenames=char('GDCh5');days=['060518'];
% mousenames=char('GDCh16');days=['090119'];
% mousenames=char('GDCh17');days=['090119'];
% mousenames=char('GDCh18');days=['080819'];
% mousenames=char('GDCh19');days=['110819'];
% mousenames=char('GDCh20');days=['080819'];
% mousenames=char('GDCh21');days=['110819'];


daynum=size(days,1);


epochl=4;
fs=256;
% calculate number of samples within an exactly 24 h long recording (not sure why this is done so complicatedly)
fsh24=fs*epochl*maxep;
output=zeros(fsh24,1);

for mouse=1:1
    mousename=mousenames(mouse,:);
    for aa=1:daynum
        
        fnout2=[mousename,'-',event,'-',days(aa,:)];
        load ([pathsig,fnout2,'.mat'],'resampled_sig')
        
        if length(resampled_sig)>fsh24
            signal=resampled_sig(1:fsh24);
        else
            signal=zeros(1,fsh24);
            signal(1:length(resampled_sig))=resampled_sig;
        end
        
        sigstim=zeros(1,length(signal));
        sigstim(find(signal>1))=1;
        sigdur=length(sigstim*2);
        
%         figure
%         subplot(4,1,aa)
%         plot(sigstim);
%          hold on

        mindur=fs;
        ba=fs;
        [startend epidur]=StimDetect(sigstim,sigdur,mindur,ba);
        
        figure;plot(epidur)
        
        figure
        scatter(startend(:,1),sigstim(startend(:,1)),'dr','filled');
        hold on
        scatter(startend(:,2),sigstim(startend(:,2)),'sm','filled');
        xlim([0 length(signal)])
%         pause
        
        fn1=[mousename,'_',days(aa,:),'_stimtest'];
        eval(['save ',pathstim,fn1,'.mat sigstim startend -mat']);
        
    end
end