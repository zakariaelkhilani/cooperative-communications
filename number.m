I=10000;
N=12;
R=2;
gss=['k-s' 'b-o' 'r-''g-*'];
%semilogy(vector_snrm,prob_outage_best,gss(:,:));hold on
%Best Relay

for snrm=1:30
channel_outage=0;
for aux=1:I

r1=sqrt(0.5)*(randn(N,1)+j*randn(N,1)); %Source-Relay Channels
h1=sqrt(0.5)*(randn(N,1)+j*randn(N,1)); %Relay Channels-Destination
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
prob_outage_best(snrm)=channel_outage/I;
vector_snrm(snrm)=snrm;
end
semilogy(vector_snrm,prob_outage_best,'-ko');
title('\fontname{Arial}outage probability using different number of relays ','FontSize',12)
xlabel('SNR dB')
ylabel('outage probability')
hold on
grid on
legend('2 relays','5 relays','8 relays','12 relays')