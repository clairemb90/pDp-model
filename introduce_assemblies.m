%strength of connections within assemblies, to adjust
pE=25; %number of connections one E assembly neurons receives from other E assembly neurons 
pI=13; %number of connections one E assembly neurons receives from I assembly neurons 
pEI=25; %number of connections one I assembly neurons receives from E assembly neurons (for tuned I, set to 0)

tit=0;

for matlist=1:20
    
    num_pat=15; % number of memories, can be changed
    num_exc=100; % E assembly size, can be changed
    num_inh=25; % I assembly size, can be changed
    
    num_ob=150;
    load(strcat('connec_A_',num2str(matlist),'.mat')) %name of the random connectivity matrix
    WE_rand=w_ee;
    WEI_rand=w_ei;
    WIE_rand=w_ie;
    load('odorseqD21_poisson.mat','OBinp') %OB input
    EnsE=[];
    Ee=zeros(num_pat,num_exc);
    pat_inc=1;
    TOT=[1:4000];
    a=1;
    for od=1:num_pat
        
        act_mc=OBinp(od,1:num_ob,1); %activated mitral cells which define odor od (150 cells with the highest firing rate)
        OB=[];
        [E_temp_l1 id_l1]=sort(sum(w_oe(act_mc(1:num_ob),:),1));
        num_conn=E_temp_l1(4000-num_exc+1);
        id_l1=id_l1(E_temp_l1>=num_conn);
        id_l1=id_l1(randperm(length(id_l1)));
        
        Ee(od,:)=id_l1(1:num_exc);
        Eout(od,:)=TOT(~ismember(TOT,Ee(od,:)));
        EnsE=[EnsE;Ee(od,:)'];
        
        for ke=1:length(Ee)
            
            ens_r=find(w_ee(Ee(od,:),Ee(od,ke))>0);
            noens_r=find(w_ee(Eout(od,:),Ee(od,ke))==1);
            if length(ens_r)<pE
                ind_b=[1:num_exc];
                n_l1=randperm(num_exc,round(pE-length(ens_r)));
                w_ee(Ee(od,ind_b(n_l1)),Ee(od,ke))=w_ee(Ee(od,ind_b(n_l1)),Ee(od,ke))+1;
                
                Eeout=Eout(od,:);
                n_l2=randperm(length(noens_r),round(pE-length(ens_r)));
                w_ee(Eeout(noens_r(n_l2)),Ee(od,ke))=0;
            else
                a=a+1;
            end
        end
        pat_inc=pat_inc+1;
        
    end
    clear Eout
    EnsE=unique(EnsE);
    NonEnsE=TOT(~ismember(TOT,EnsE));
    
    EnsI=[];
    Ii=zeros(num_pat,num_inh);
    
    pat_inc=1;
    
    for od=1:num_pat
        
        a=1;
        b=1;
        
        [I_pre id_i]=sort(sum(w_ei(Ee(od,:),:),1));
        
        num_conn_I=I_pre(length(id_i)-num_inh+1);
        id_i=id_i(I_pre>=num_conn_I);
        id_i=id_i(randperm(length(id_i)));
        Ii(od,:)=id_i(1:num_inh);
        TOTI=[1:1000]; %FB
        Iout(od,:)=TOTI(~ismember(TOTI,Ii(od,:)));
        
        EnsI=[EnsI;Ii(od,:)'];
        
        for ke=1:length(Ee)
            ens_rI=find(w_ie(Ii(od,:),Ee(od,ke))>0);
            noens_rI=find(w_ie(Iout(od,:),Ee(od,ke))==1);
            if length(ens_rI)<pI
                ind_bI=[1:num_inh];
                n_l3=randperm(num_inh,round(pI-length(ens_rI)));
                w_ie(Ii(od,ind_bI(n_l3)),Ee(od,ke))=w_ie(Ii(od,ind_bI(n_l3)),Ee(od,ke))+1;
                
                Iiout=Iout(od,:);
                n_l4=randperm(length(noens_rI),round(pI-length(ens_rI)));
                w_ie(Iiout(noens_rI(n_l4)),Ee(od,ke))=0;
            else
                b=b+1;
            end
        end
        
        TOTE=[1:4000];
        Eout(od,:)=TOTE(~ismember(TOTE,Ee(od,:)));
        
        if pEI~=0
            for ki=1:length(Ii)
                
                ens_rI=find(w_ei(Ee(od,:),Ii(od,ki))>0);
                noens_rI=find(w_ei(Eout(od,:),Ii(od,ki))==1);
                
                if length(ens_rI)<pEI
                    n_l3=randperm(length(Ee(od,:)),round(pEI-length(ens_rI)));
                    ind_bE=[1:length(Ee(od,:))];
                    w_ei(Ee(od,ind_bE(n_l3)),Ii(od,ki))=w_ei(Ee(od,ind_bE(n_l3)),Ii(od,ki))+1;
                    
                    Iiout=Eout(od,:);
                    n_l4=randperm(length(noens_rI),round(pEI-length(ens_rI)));
                    w_ei(Iiout(noens_rI(n_l4)),Ii(od,ki))=0;
                end
                
            end
            pat_inc=pat_inc+1;
        end
        
    alpha_od(od)=mean(sum(w_ee(Ee(od,:),Ee(od,:))));
    beta_od(od)=mean(sum(w_ie(Ii(od,:),Ee(od,:))));
    gamma_od(od)=mean(sum(w_ei(Ee(od,:),Ii(od,:))));
    alpha_ref(od)=mean(sum(WE_rand(Ee(od,:),Ee(od,:))));
    beta_ref(od)=mean(sum(WIE_rand(Ii(od,:),Ee(od,:))));
    gamma_ref(od)=mean(sum(WEI_rand(Ee(od,:),Ii(od,:))));
    
    end

    alpha=mean(alpha_od./alpha_ref);
    beta=mean(beta_od./beta_ref);
    gamma=mean(gamma_od./gamma_ref);
    
    save(strcat('connec_A_',num2str(matlist),'_tEI_noA.mat'),'pEI','pI','pE','w_oe','w_oi','w_ie','w_ei','w_ii','w_ee','EnsE','Ee','Ii','EnsI','num_pat','alpha','beta','gamma')
    clearvars -except matlist pE pI pEI tit
    
end