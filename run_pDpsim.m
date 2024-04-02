tit=1;

twoorthree=input('network architecture: one inhibitory population (1) or two inhibitory populations, one feedback and one feedworward (2) ');
randortuned=input('network type: random (1), scaled (2), tuned E+I (3) or tuned I (4) ');
analys=input('analysis: yes (1) or no (2) ');

if twoorthree==1
    pagexls=1;
    if randortuned==1
        ext='';
    elseif randortuned==2
        ext='';
        ext1='_tI';
        fact=1.4;
    elseif randortuned==3
        ext='_tEI';
    elseif randortuned==4
        ext='_tI';
    end
    ods=input('odor set: learned (1), novel (2), subspace fig.4 (3) or subspace fig. 6 (4) ');
    if ods==1
        num_odors=10;
        load('odorset_10learned.mat')
    elseif ods==2
        num_odors=10;
        load('odorset_10novel.mat')
    elseif ods==3
        load('odorspace_fig4.mat')
    elseif ods==4
        load('odorspace_fig6.mat')
    else
        disp('error')
    end
    PIN=0;
else
    pagexls=2;
    load('odorset_20odors.mat');num_odors=20;
    PIN=input('inhibiting I neurons (vPIN): no (0), FFI (1) or FBI (7)');
    if randortuned==1
        ext='_r';
    elseif randortuned==2
        ext='_r';
        ext1='_s';
        fact=disp('scaling factor I to E strength?');
    elseif randortuned==4
        ext='_s';
    elseif randortuned==3
        disp('network not available')
    end
end

[param,connec,tot]=xlsread('parameters.xlsx',pagexls);
if twoorthree==2 && randortuned==1
    tot_mat=10;
else
    tot_mat=length(param);
end

r_ob=int8(r_olfbs);clear r_olfbs

for mat=1:tot_mat %can be changed if simulation of only subsets of networks    
   
    Ninh=1000;
    if twoorthree==2
        Istim=ones(length(r_ob),Ninh);
        PA1=0;
        if PIN>0
            [inh_neur]=xlsread('parameters.xlsx',3); %neurons identity that are silenced can be changed: these are the ones used in the paper
            for oo=1:num_odors
                times=[500/dt+3000/dt*(oo-1);3000/dt+3000/dt*(oo-1)];
                Istim(times(1):times(2),[inh_neur(oo,PIN):inh_neur(oo,PIN+1),inh_neur(mat,PIN+2):inh_neur(mat,PIN+3)])=PA1;
                clear times
            end
        end
    end
    
    struct=1;
    
    if randortuned==2
        load(strcat(connec{mat},ext1)),load(strcat(connec{mat},ext),'w_ie');
        w_ie=fact*w_ie;
    else
        load(strcat(connec{mat},ext))
    end
    
    Nexc=4000;
    Ninh=length(w_ii);
    Nob=1500;

    y=param(mat,:);
    
    dt=0.1;
    
    %neuronal parameters
    tau_exc=30;
    tau_inh=10;
    E_e=0; 
    E_i=-70; 
    g_rest_e=1.35;
    g_rest_i=0.9;
    tau_e=85; 
    tau_i=50;
    tau_ref_e=8/dt;
    tau_ref_i=8/dt;
    tau_w=40;
    ada=1;
   
    E_th=-38; 
    E_th_i=-45;
    
    if twoorthree==1
        E_rest=-60; %
        E_rest_ii=-65;
        beka=10;
    elseif twoorthree==2
        E_rest=-60;
        E_rest_ii=-60;
        beka=5;
    end
    
    %activity neuron
    if ods>2
        ds=0.01;
    else
        ds=0.1;
    end
    downs=floor(length(r_ob)*ds);
    perc=1;
    NumNe=[1:Nexc/perc];
    NumNi=[1:Ninh/perc];
    V_id=E_rest_ii.*ones(downs,Ninh/perc);
    V_ed=E_rest.*ones(downs,Nexc/perc);
    g_oed=zeros(downs,Nexc/perc);
    g_eed=zeros(downs,Nexc/perc);
    g_ied=zeros(downs,Nexc/perc);
    g_ief=zeros(downs,Nexc/perc);
    g_ieb=zeros(downs,Nexc/perc);
   
    %% solve
    t=1;
    ts=1;
   
    tsc=1;
    tsci=1;
    tsci_b=1;
    VI=E_rest_ii.*ones(1,Ninh);
    VE=E_rest.*ones(1,Nexc);
    GEE=zeros(1,Nexc);
    GIE=zeros(1,Nexc);
    GIEB=zeros(1,Nexc);
    GIEF=zeros(1,Nexc);
    GEI=zeros(1,Ninh);
    GII=zeros(1,Ninh);
    GOE=zeros(1,Nexc);
    GOI=zeros(1,Ninh);
    
    refrac_E=single(zeros(1,Nexc));
    refrac_I=single(zeros(1,Ninh));
    WE=zeros(1,Nexc);
  
    spikecount_E=zeros(1000000,2);
    spikecount_I=zeros(3000000,2);
    
    while t<=length(r_ob)
        
        V_e=VE;
        V_i=VI;
        g_ee=GEE;
        g_ibe=GIEB;
        g_ife=GIEF;
        if twoorthree==1
            g_ie=GIE;
        else
            g_ie=GIEB+GIEF;
        end
        g_ei=GEI;
        g_ii=GII;
        g_oe=GOE;
        g_oi=GOI;
        w_e=WE;
        
        VE=V_e+dt/tau_e.*((E_rest-V_e)+((g_oe+g_ee).*(E_e-V_e)+g_ie.*(E_i-V_e)-w_e)/g_rest_e);
        
        VI=V_i+dt/tau_i.*((E_rest_ii-V_i)+((g_oi+g_ei).*(E_e-V_i)+g_ii.*(E_i-V_i))/g_rest_i);
        
        WE=w_e+dt/tau_w*(ada*(V_e-E_rest)-w_e);
        
        refre=find(refrac_E(1,:)>0);
        refri=find(refrac_I(1,:)>0);
        VE(refre)=E_rest;
        VI(refri)=E_rest_ii;
        neufire=find(VE>E_th);
        neufiri=find(VI>E_th_i);
        
        spikecount_E(tsc:tsc+length(neufire)-1,1)=t;
        spikecount_E(tsc:tsc+length(neufire)-1,2)=neufire;
        tsc=tsc+length(neufire);
        spikecount_I(tsci:tsci+length(neufiri)-1,1)=t;
        spikecount_I(tsci:tsci+length(neufiri)-1,2)=neufiri;
        tsci=tsci+length(neufiri);
        obfire=find(r_ob(t,:)>0);
        refrac_E=refrac_E-1;
        refrac_I=refrac_I-1;
        refrac_E(neufire)=tau_ref_e;
        refrac_I(neufiri)=tau_ref_i;
        T=2;
        WE(1,neufire)=WE(1,neufire)+beka;
       
        GEE=g_ee.*exp(-dt/tau_exc)+y(1).*sum(w_ee(neufire,:),1);
        GEI=g_ei.*exp(-dt/tau_exc)+y(4).*sum(w_ei(neufire,:),1);
        if twoorthree==1
            GIE=g_ie.*exp(-dt/tau_inh)+y(2).*sum(w_ie(neufiri,:),1);
            GII=g_ii.*exp(-dt/tau_inh)+y(6).*sum(w_ii(neufiri,:),1);
        else
            GIEB=g_ibe.*exp(-dt/tau_inh);
            GIEF=g_ife.*exp(-dt/tau_inh);
            for in=1:length(neufiri)
                if neufiri(in)<Ninh/2+1
                    GIEB=GIEB+y(7).*squeeze(Istim(t,neufiri(in)))*w_ie(neufiri(in),:);
                else
                    GIEF=GIEF+y(2).*squeeze(Istim(t,neufiri(in)))*w_ie(neufiri(in),:);
                end
            end
            GIE=GIEB+GIEF;
            GII=g_ii.*exp(-dt/tau_inh)+y(6).*sum(squeeze(Istim(t,neufiri))*w_ii(neufiri,:),1);
        end
        
        GOE=g_oe.*exp(-dt/tau_exc)+y(3).*sum(w_oe(obfire,:),1);
        GOI=g_oi.*exp(-dt/tau_exc)+y(5).*sum(w_oi(obfire,:),1);
        if mod(t-1,1/ds)==0 || t==1
            V_id(ts,:)=VI(1,NumNi);
            V_ed(ts,:)=VE(1,NumNe);
        
            g_oed(ts,:)=GOE(1,NumNe);
            g_eed(ts,:)=GEE(1,NumNe);
            g_ied(ts,:)=GIE(1,NumNe);
            if twoorthree==2
                g_ieb(ts,:)=GIEB(1,NumNe);
                g_ief(ts,:)=GIEF(1,NumNe);
            end
            ts=ts+1;
        end
        t=t+1;
        clear refre refri neufire neufiri obfire
    end
    clear refrac_E refrac_I 
    spikecount_E=spikecount_E(1:tsc-1,:);
    spikecount_I=spikecount_I(1:tsci-1,:);
    g_ed=g_oed+g_eed;
    
    %% ANALYSIS
    ACT{tit}=spikecount_E;
    ACT_I{tit}=spikecount_I;
    
    if analys==1
        observables
        co_tuning
    end
    tit=tit+1
    clearvars -except mat twoorthree tit randortuned analys pagexls ext ext1 fact ods num_odors r_ob tot_mat PIN ACT ACT_I param connec tot
    %save()
end

