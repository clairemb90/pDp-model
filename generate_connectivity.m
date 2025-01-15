function [w_oe,w_oi,w_ee,w_ii,w_ei,w_ie]=generate_connectivity(p_oe,p_oi,p_ee,p_ei,p_ief,p_ieb,p_ii_wb,p_ii_wf,num,typ)
%inputs:

% -typ: 1 for 1 type of inhibitory population, 2 for feedback+feedforward
%       inhibitory populations
% -p_XY:probabilities of connections from presynaptic neuron X to
%       postsynaptic neuron Y p_XY. If 1 type of inhibitory population,
%       p_ii is given by p_ii_wb and p_ie by p_ieb
% -num: file saved as connec_rand_"number defined in num".mat

Ninh=1000;
Nexc=4000;
Nob=1500;
Ifb=500; %half feedback, half feedforward

if typ==1
    p_ii=p_ii_wb;
    p_ie=p_ieb;

    for ne=1:Nexc
        w_oe(:,ne)=randperm(Nob,Nob)/Nob;
        w_ee(:,ne)=randperm(Nexc,Nexc)/Nexc;
        w_ie(:,ne)=randperm(Ninh,Ninh)/Ninh;
    end
    for ni=1:Ninh
        w_oi(:,ni)=randperm(Nob,Nob)/Nob;
        w_ei(:,ni)=randperm(Nexc,Nexc)/Nexc;
        w_ii(:,ni)=randperm(Ninh,Ninh)/Ninh;
    end
else
    for ne=1:Nexc
        w_oe(:,ne)=randperm(Nob,Nob)/Nob;
        w_ee(:,ne)=randperm(Nexc,Nexc)/Nexc;
        w_ieb(1:Ifb,ne)=randperm(Ifb,Ifb)./(Ifb);
        w_ief(Ifb+1:Ninh,ne)=randperm((Ninh-Ifb),(Ninh-Ifb))./((Ninh-Ifb));
    end
    for ni=1:Ifb
        w_ei(:,ni)=randperm(Nexc,Nexc)/4000;
        w_ii_wb(1:Ifb,ni)=randperm(Ifb,Ifb)./(Ifb);

    end
    for ni=Ifb+1:Ninh
        w_oi(:,ni)=randperm(Nob,Nob)/Nob;
        w_ii_wf(Ifb+1:Ninh,ni)=randperm((Ninh-Ifb),(Ninh-Ifb))./((Ninh-Ifb));
    end
end

w_oe(w_oe>=1-p_oe)=1;
w_oe(w_oe<1)=0;

w_ee(w_ee>=1-p_ee)=1;
w_ee(w_ee<1)=0;

w_oi(w_oi>=1-p_oi)=1;
w_oi(w_oi<1)=0;

w_ei(w_ei>=1-p_ei)=1;
w_ei(w_ei<1)=0;

if typ==1
    w_ie(w_ie>=1-p_ie)=1;
    w_ie(w_ie<1)=0;

    w_ii(w_ii>=1-p_ii)=1;
    w_ii(w_ii<1)=0;
else
    w_oi(:,1:Ifb)=0;
    w_ei(:,Ifb+1:Ninh)=0;

    w_ii_wb(w_ii_wb>=1-p_ii_wb)=1;
    w_ii_wb(w_ii_wb<1)=0;

    w_ii_wf(w_ii_wf>=1-p_ii_wf)=1;
    w_ii_wf(w_ii_wf<1)=0;
    w_ii(1:Ifb,1:Ifb)=w_ii_wb(1:Ifb,1:Ifb);

    w_ii(Ifb+1:Ninh,Ifb+1:Ninh)=w_ii_wf(Ifb+1:Ninh,Ifb+1:Ninh);

    w_ii(1:Ifb,Ifb+1:Ninh)=0;
    w_ii(Ifb+1:Ninh,1:Ifb)=0;

    w_ieb(w_ieb>=1-p_ieb)=1;
    w_ieb(w_ieb<1)=0;

    w_ief(w_ief>=1-p_ief)=1;
    w_ief(w_ief<1)=0;

    w_ie(1:Ifb,:)=w_ieb(1:Ifb,:);
    w_ie(Ifb+1:Ninh,:)=w_ief(Ifb+1:Ninh,:);

    %
end

TOTO=[1:Nob];
TOTE=[1:Nexc];

%lower the variation in the number of output connections across neurons
ee_p=[find(sum(w_ee,2)>p_ee*Nexc*1.03);  find(sum(w_ee,2)<p_ee*Nexc*0.97)];
s_ep=length(ee_p);
for p=1:s_ep
    ee_n=find(w_ee(ee_p(p),:)==1)';
    ee_r=TOTE(~ismember(TOTE,ee_n));
    if length(ee_n)>p_ee*Nexc*1.03
        w_ee(ee_p(p),ee_n(randperm(length(ee_n),round(length(ee_n)-p_ee*Nexc*1.03))))=0;
    elseif length(ee_n)<p_ee*Nexc*0.97
        w_ee(ee_p(p),ee_r(randperm(length(ee_r),round(p_ee*Nexc*0.97-length(ee_n)))))=1;
    end
end

if typ==1
    TOTI=[1:Ninh];
    ie_p=[find(sum(w_ie,2)>p_ie*Nexc*1.05);  find(sum(w_ie,2)<p_ie*Nexc*0.95)];
    s_iep=length(ie_p);
    for p=1:s_iep
        iee_n=find(w_ie(ie_p(p),:)==1)';
        iee_r=TOTE(~ismember(TOTE,iee_n));
        if length(iee_n)>p_ie*Nexc*1.03
            w_ie(ie_p(p),iee_n(randperm(length(iee_n),round(length(iee_n)-p_ie*Nexc*1.03))))=0;
        elseif length(iee_n)<p_ie*Nexc*0.97
            w_ie(ie_p(p),iee_r(randperm(length(iee_r),round(p_ie*Nexc*0.97-length(iee_n)))))=1;
        end
    end

    ei_p=[find(sum(w_ei,2)>p_ei*Ninh*1.05);  find(sum(w_ei,2)<p_ei*Ninh*0.95)];
    s_eip=length(ei_p);
    for p=1:s_eip
        eie_n=find(w_ei(ei_p(p),:)==1)';
        eie_r=TOTI(~ismember(TOTI,eie_n));
        if length(eie_n)>p_ei*Ninh*1.03
            w_ei(ei_p(p),eie_n(randperm(length(eie_n),round(length(eie_n)-p_ei*Ninh*1.03))))=0;
        elseif length(eie_n)<p_ei*Ninh*0.97
            w_ei(ei_p(p),eie_r(randperm(length(eie_r),round(p_ei*Ninh-length(eie_n)))))=1;
        end
    end

    ii_p=[find(sum(w_ii,2)>p_ii*Ninh*1.03);  find(sum(w_ii,2)<p_ii*Ninh*0.97)];
    s_ip=length(ii_p);
    for p=1:s_ip
        ii_n=find(w_ii(ii_p(p),:)==1)';
        ii_r=TOTI(~ismember(TOTI,ii_n));
        if length(ii_n)>p_ii*Ninh*1.03
            w_ii(ii_p(p),ii_n(randperm(length(ii_n),round(length(ii_n)-p_ii*Ninh*1.03))))=0;
        elseif length(ii_n)<p_ii*Ninh*0.97
            w_ii(ii_p(p),ii_r(randperm(length(ii_r),round(p_ii*Ninh*0.97-length(ii_n)))))=1;
        end
    end
else
    TOTI=[1:500];
    WEIb=w_ei(:,1:500);
    p_ei=mean(sum(WEIb,2))/500;

    ei_p=[find(sum(WEIb,2)>p_ei*Ninh/2*1.04);  find(sum(WEIb,2)<p_ei*Ninh/2*0.96)];
    s_eip=length(ei_p);
    for p=1:s_eip
        eie_n=find(WEIb(ei_p(p),:)==1)';
        eie_r=TOTI(~ismember(TOTI,eie_n));
        if length(eie_n)>p_ei*Ninh/2*1.04
            WEIb(ei_p(p),eie_n(randperm(length(eie_n),round(length(eie_n)-p_ei*Ninh/2*1.04))))=0;
        elseif length(eie_n)<p_ei*Ninh/2*0.96
            WEIb(ei_p(p),eie_r(randperm(length(eie_r),round(p_ei*Ninh/2*0.96-length(eie_n)))))=1;
        end
    end
    w_ei(:,1:500)=WEIb;
    %
    WIEb=w_ie(1:500,:);
    ie_p=[find(sum(WIEb,2)>p_ieb*Nexc*1.05);  find(sum(WIEb,2)<p_ieb*Nexc*0.95)];
    s_iep=length(ie_p);
    for p=1:s_iep
        iee_n=find(WIEb(ie_p(p),:)==1)';
        iee_r=TOTE(~ismember(TOTE,iee_n));
        if length(iee_n)>p_ieb*Nexc*1.03
            WIEb(ie_p(p),iee_n(randperm(length(iee_n),round(length(iee_n)-p_ieb*Nexc*1.03))))=0;
        elseif length(iee_n)<p_ieb*Nexc*0.97
            WIEb(ie_p(p),iee_r(randperm(length(iee_r),round(p_ieb*Nexc*0.97-length(iee_n)))))=1;
        end
    end

    WIEf=w_ie(501:1000,:);
    ie_p=[find(sum(WIEf,2)>p_ief*Nexc*1.05);  find(sum(WIEf,2)<p_ief*Nexc*0.95)];
    s_iep=length(ie_p);
    for p=1:s_iep
        iee_n=find(WIEf(ie_p(p),:)==1)';
        iee_r=TOTE(~ismember(TOTE,iee_n));
        if length(iee_n)>p_ief*Nexc*1.03
            WIEf(ie_p(p),iee_n(randperm(length(iee_n),round(length(iee_n)-p_ief*Nexc*1.03))))=0;
        elseif length(iee_n)<p_ief*Nexc*0.97
            WIEf(ie_p(p),iee_r(randperm(length(iee_r),round(p_ief*Nexc*0.97-length(iee_n)))))=1;
        end
    end

    w_ie(1:500,:)=WIEb;
    w_ie(501:1000,:)=WIEf;

    TOTI=[1:Ninh/2];
    WIIb=w_ii(1:500,1:500);
    p_ii=p_ii_wb;
    ii_p=[find(sum(WIIb,2)>p_ii*Ninh/2*1.02);  find(sum(WIIb,2)<p_ii*Ninh/2*0.98)];
    s_ip=length(ii_p);
    for p=1:s_ip
        ii_n=find(WIIb(ii_p(p),:)==1)';
        ii_r=TOTI(~ismember(TOTI,ii_n));
        if length(ii_n)>p_ii*Ninh/2*1.02
            WIIb(ii_p(p),ii_n(randperm(length(ii_n),round(length(ii_n)-p_ii*Ninh/2*1.02))))=0;
        elseif length(ii_n)<p_ii*Ninh/2*0.98
            WIIb(ii_p(p),ii_r(randperm(length(ii_r),round(p_ii*Ninh/2*0.98-length(ii_n)))))=1;
        end
    end

    WIIf=w_ii(501:1000,501:1000);
    p_ii=p_ii_wf;
    ii_p=[find(sum(WIIf,2)>p_ii*Ninh/2*1.02);  find(sum(WIIf,2)<p_ii*Ninh/2*0.98)];
    s_ip=length(ii_p);
    for p=1:s_ip
        ii_n=find(WIIf(ii_p(p),:)==1)';
        ii_r=TOTI(~ismember(TOTI,ii_n));
        if length(ii_n)>p_ii*Ninh/2*1.02
            WIIf(ii_p(p),ii_n(randperm(length(ii_n),round(length(ii_n)-p_ii*Ninh/2*1.02))))=0;
        elseif length(ii_n)<p_ii*Ninh/2*0.98
            WIIf(ii_p(p),ii_r(randperm(length(ii_r),round(p_ii*Ninh/2*0.98-length(ii_n)))))=1;
        end
    end
    w_ii(1:500,1:500)=WIIb;
    w_ii(501:1000,501:1000)=WIIf;

end

save(strcat('connec_E_',num2str(num),'.mat'),'w_oe','w_oi','w_ee','w_ii','w_ei','w_ie','p_ee','Nexc','Ninh','Nob');
%
figure,
subplot(3,3,[1 2 4 5])
imagesc(w_ee)
subplot(3,3,[3 6])
plot(sum(w_ee,1),fliplr(1:Nexc))
subplot(3,3,[7 8])
plot(sum(w_ee,2))


