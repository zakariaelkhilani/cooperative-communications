clear all
%close all
I=50000;
N=10;
R=2;
% %Alamouti
for snrm=1:30
 channel_outage=0;
for aux=1:I
r1=sqrt(0.5)*(randn(N,1)+j*randn(N,1)); 
h1=sqrt(0.5)*(randn(N,1)+j*randn(N,1)); %relay to distination 
snrr1=(abs(r1).^2)*snrm; %SNR 
vect_outage=snrr1>=(2^(2*R)-1);
if(sum(vect_outage)>=2)
 vect_sec_link=sort(abs(h1).^2.*vect_outage,'descend');
 snrr2=(snrm/2)*(sum(vect_sec_link(1:2)));
 channel_outage=channel_outage+(snrr2<(2^(2*R)-1));
else
 channel_outage=channel_outage+1;
end
end
vector_snrm(snrm)=snrm;
prob_outage_Alamouti(snrm)=channel_outage/I;
end

%Best Relay
for snrm=1:30
 channel_outage=0;
for aux=1:I

r1=sqrt(0.5)*(randn(N,1)+j*randn(N,1)); %Canales Fuente-Relay
h1=sqrt(0.5)*(randn(N,1)+j*randn(N,1)); %Canales Relay-Destino
snrr1=(abs(r1).^2)*snrm;
vect_outage=snrr1>=(2^(2*R)-1);
if(sum(vect_outage)>0)
 vect_sec_link=sort(abs(h1).^2.*vect_outage,'descend');
 snrr2=snrm*vect_sec_link(1);
 channel_outage=channel_outage+(snrr2<(2^(2*R)-1));
else
 channel_outage=channel_outage+1;
end
end
vector_snrm(snrm)=snrm;
prob_outage_Best(snrm)=channel_outage/I;
end
%Direct transmission
for snrm=1:30
 channel_outage=0;
for aux=1:I

r1=sqrt(0.5)*(randn(1,1)+j*randn(1,1)); %Source-Destination Channels
 snrr2=snrm*(abs(r1).^2);
 channel_outage=channel_outage+(snrr2<(2^(2*R)-1));
end
prob_outage_Direct(snrm)=channel_outage/I;
vector_snrm(snrm)=snrm;
end


semilogy(vector_snrm,prob_outage_Direct,'g-o');
hold on
grid on
semilogy(vector_snrm,prob_outage_Best,'r-*');
hold on
grid on
semilogy(vector_snrm,prob_outage_Alamouti,'b->');

%title('\fontname{Arial}Outage probability when distances are taken into account ','FontSize',12)
xlabel('SNR dB')
 ylabel('Outage Probability');
 hold on
 axis([0 30 10^-4 1])
 legend('Direct link','Best relay','Alamouti','Alamouti dis')