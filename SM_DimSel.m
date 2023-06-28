%% This code is for Optimum subspace dimension selection
% This code validates that the maximum classification accuracy is always
% between dimension 1 and the dimension of maximum entry in overlapping
% matrix
% Considered only one session of data of one subject, which can be 
% generalized to any combination of sessions. Used 10-fold cross validation 
% to validate classification performance.

clear all;
close all;
clc;

fs=250;
addpath([pwd '/../dataset/BCI4_2B/']);

% BCI competition 4 Dataset 2B
[s_sess3,h_sess3]=sload('B0403T.gdf'); % B0403T dataset preparation
Trials_All_Ind_sess3=find(h_sess3.EVENT.TYP==768);
Reject_Ind_sess3=find(h_sess3.ArtifactSelection==1);

Trials_All_Ind_sess3(Reject_Ind_sess3)=[];
h_sess3.Classlabel(Reject_Ind_sess3)=[];

Trials_Ind_sess3=Trials_All_Ind_sess3;
Z_sess3= h_sess3.Classlabel';
trial=1;
for i=1:length(Trials_Ind_sess3)
    st_sess3=h_sess3.EVENT.POS(Trials_Ind_sess3(i));
    dur_sess3=h_sess3.EVENT.DUR(Trials_Ind_sess3(i));
    Y_temp_sess3(:,:)=s_sess3(st_sess3:st_sess3+1999,1:3);
    if(find(isnan(Y_temp_sess3(:,:)))~=0) 
       i
       Z_sess3(trial)=[];
    else
        Y_sess3(:,:,trial)=Y_temp_sess3(:,:)';
        trial=trial+1;
    end
end

st=3.5*fs % Considering only 3.5 to 5.5 seconds duration of each trial
ed=5.5*fs
Ch=3;

ZTn=[Z_sess3];
YTn=cat(3,Y_sess3);

indL=find(ZTn==1);
indR=find(ZTn==2);
NL=length(indL);
NR=length(indR);
L=min([NL,NR]);

for d=1:2*Ch-1
    clear BPL; clear BPL_Temp;  clear BPR; clear BPR_Temp; 
    clear Ur; clear Ul; 
    % Bandpower features calculation
    for trial=1:L
        Ytemp(:,:)=YTn(:,st:ed,indL(trial));
        BP_Temp(:,:)=[bandpower(Ytemp',fs,[8 12]);bandpower(Ytemp',fs,[12 30])];
        BPL_Temp=[];
        for c=1:Ch
            BPL_Temp=[BPL_Temp;BP_Temp(:,c)];
        end
        BPL(:,trial)=BPL_Temp;
    end
    
    for trial=1:L
        Ytemp(:,:)=YTn(:,st:ed,indR(trial));
        BP_Temp(:,:)=[bandpower(Ytemp',fs,[8 12]);bandpower(Ytemp',fs,[12 30])];
        BPR_Temp=[];
        for c=1:Ch
            BPR_Temp=[BPR_Temp;BP_Temp(:,c)];
        end
        BPR(:,trial)=BPR_Temp;
    end
    c = cvpartition(L,'KFold',10);
    clear O;
    tDs=d;
    for fold=1:c.NumTestSets
        clear test_ind; clear train_ind; clear Rr; clear Rl; clear Ur; clear Er; clear Ul; clear El;
        test_ind=find(test(c,fold)==1); 
        train_ind=find(training(c,fold)==1);

        % Subspace learning of train data
        Rr=cov(BPR(:,train_ind)');  Rl=cov(BPL(:,train_ind)'); 
        [Ur Er ~]=svd(Rr); [Ul El ~]=svd(Rl); 
   
        % Computing pairwise ovelapping matrix
        for pi=1:6
            for pj=1:6
                O(pi,pj,fold)=subspace(Ul(:,1:pi),Ur(:,1:pj));
            end
        end
        clear RrTt; clear ErTt; clear UrTtl; clear ErTtl; clear UrTt; clear Zpredict; clear ZTrue; clear ZpredictL; clear ZpredictR;
        for i=1:length(test_ind)
            RrTt=BPR(:,test_ind(i))*BPR(:,test_ind(i))'; 
            [UrTt ErTt ~]=svd(RrTt); 
            [UrTtl ErTtl ~]=svd(UrTt(:,1:d)'*Ul(:,1:d));
            cos_theta=diag(ErTtl);
            SrTtl=sum(cos_theta(1:tDs).^2)/tDs;
    
            [UrTtr ErTtr ~]=svd(UrTt(:,1:d)'*Ur(:,1:d));
            cos_theta=diag(ErTtr);
            SrTtr=sum(cos_theta(1:tDs).^2)/tDs;
            [VV ZpredictR(i)]=max([SrTtl SrTtr]);
           
            clear RlTt; clear UlTt; clear ElTt; clear UlTtl; clear ElTtl; clear UlTtr; clear ElTtr;
            RlTt=BPL(:,test_ind(i))*BPL(:,test_ind(i))';
            [UlTt ElTt ~]=svd(RlTt); 
            [UlTtl ElTtl ~]=svd(UlTt(:,1:d)'*Ul(:,1:d));
            cos_theta=diag(ElTtl);
            SlTtl=sum(cos_theta(1:tDs).^2)/tDs;
    
            [UlTtr ElTtr ~]=svd(UlTt(:,1:d)'*Ur(:,1:d));
            cos_theta=diag(ElTtr);
            SlTtr=sum(cos_theta(1:tDs).^2)/tDs;
            [VV ZpredictL(i)]=max([SlTtl SlTtr]);
       end
       ZTrue=[zeros(1,length(test_ind))+1 zeros(1,length(test_ind))+2];
       Zpredict=[ZpredictL ZpredictR];
       [Cf order]=confusionmat(ZTrue,Zpredict);
       Acc(fold)=(Cf(1,1)+Cf(2,2))/length(ZTrue);    
   end
   Acc_d(d)=mean(Acc);
end
contour(mean(O,3),'ShowText','on','Linewidth',3);
title('Contour plot of Overlapping matrix')