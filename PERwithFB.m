clear all;
%close all;
I=100000;
N=10;%number of relays 
R=2;%rate 
d=100;
%Calculate the distancies of each destination
distances=d*rand(N,1);
distancey=(d/2)*rand(N,1);
distance=sqrt(distances.^2+distancey.^2);
distances2=d-distances;
distance2=sqrt(distances2.^2+distancey.^2);
fc=2.5*10^9;
c=3*10^8;
landac=c/fc;
d0=1;
mu=3;
Pn=1;
Probe=0.05;
L=10;%number of symbols in each packets 
%We chose SNR best and used Alamouti
for snrm=1:30
 channel_outage_Alamouti=0;
 rate_mean_Alamouti=0;
 rateslot=R;
 snrml=10^(snrm/10);
 Nslots=0;
 variance1(snrm)=((d/d0)^(-mu))*(landac/(4*pi*d0))^2;
 Ps=snrml*Pn/variance1(snrm);
 for aux=1:I
 for k=1:N
 variance2=((distance(k)/d0)^(-mu))*(landac/(4*pi*d0))^2;
 r(k,1)=sqrt(0.5)*sqrt(variance2)*(randn(1,1)+j*randn(1,1));%rayleigh
 end
 for k=1:N
 variance3=((distance2(k)/d0)^(-mu))*(landac/(4*pi*d0))^2;
 h(k,1)=sqrt(0.5)*sqrt(variance3)*(randn(1,1)+j*randn(1,1));
 end
 Nslots=Nslots+1;
 snrr1=(abs(r).^2)*Ps/(Pn*2);
 vect_outage=snrr1>=(2^(2*R)-1);
 if(sum(vect_outage)>=2)
 %We look at the feedback error
 vect_sec_link=abs(h).^2.*vect_outage;
 [maxim,pos1]=max(vect_sec_link); %We take the maximum
 vect_sec_link2(1)=vect_sec_link(pos1); %We keep it
 vect_sec_link(pos1)=0; %Put your position to 0
 [maxim2,pos2]=max(vect_sec_link); %We take the next maximum
 vect_sec_link(pos1)=vect_sec_link2(1); %We return the maximum value in its position
 pos1=pos1-1; %To start positions at 0
 pos2=pos2-1;
 b1=dec2bin(pos1); % We pass the positions to binary
 b2=dec2bin(pos2);
  a=rand(1);
  b=rand(1);
  c=rand(1);
  dd=rand(1);
 l=0;
 l(4)=b1(length(b1));
 
 if pos1>1
 l(3)=b1(length(b1)-1);
 end
 
 if pos1>3
 l(2)=b1(length(b1)-2);
 end

 if pos1>7
 l(1)=b1(length(b1)-3);
 end

 l=l>'0';



 if a<=Probe
 l(4)=not(l(4));

 end

 if b<=Probe && pos1>1
 l(3)=not(l(3));
 end

 if c<=Probe && pos1>3
 l(2)=not(l(2));
 end

 if dd<=Probe && pos1>7
 l(1)=not(l(1));
 end


 %Convert binary to decimal

 choice1=l(4)+2*l(3)+4*l(2)+8*l(1);



 a=rand(1);
 b=rand(1);
 c=rand(1);
 dd=rand(1);

 l=0;

 l(4)=b2(length(b2));

 if pos2>1

 l(3)=b2(length(b2)-1);
 end

 if pos2>3
 l(2)=b2(length(b2)-2);
 end

 if pos2>7
 l(1)=b2(length(b2)-3);
 end

 l=l>'0';



 if a<=Probe
 l(4)=not(l(4));

 end

 if b<=Probe && pos1>1
 l(3)=not(l(3));
 end

 if c<=Probe && pos1>3
 l(2)=not(l(2));
 end

 if dd<=Probe && pos1>7
 l(1)=not(l(1));
 end

 %Convert binary to decimal
 choice2=l(4)+2*l(3)+4*l(2)+8*l(1);

 
 choice1=choice1+1;
 choice2=choice2+1;



 if N<choice1
 vect_sec_link(choice1)=0;
 end

 if N<choice2
 vect_sec_link(choice2)=0;
 end

 snrr2=(Ps/(Pn*4))*(vect_sec_link(choice1)+vect_sec_link(choice2));

 channel_outage_Alamouti=channel_outage_Alamouti+(snrr2<(2^(2*R)-1));

 if snrr2>=(2^(2*R)-1)
 rate_mean_Alamouti=rate_mean_Alamouti+rateslot;
 else
 rate_mean_Alamouti=rate_mean_Alamouti+0;
 end


 else

 channel_outage_Alamouti=channel_outage_Alamouti+1;
 rate_mean_Alamouti=rate_mean_Alamouti+0;

 end
 end
 vector_snrm(snrm)=snrm;
 prob_outage_Alamouti(snrm)=channel_outage_Alamouti/I;
 prob_error_paquet_Alamouti(snrm)=1-(1-prob_outage_Alamouti(snrm))^(L);
 rate_Alamouti(snrm)=rate_mean_Alamouti/Nslots;
end


%Choosing only 1
for snrm=1:30
 Nslots=0;
 rate_mean_Best=0;
 channel_outage_Best=0;
 snrml=10^(snrm/10);

 variance1(snrm)=((d/d0)^(-mu))*(landac/(4*pi*d0))^2;
 Ps=snrml*Pn/variance1(snrm);
 for aux=1:I
 for k=1:N

 variance2=((distance(k)/d0)^(-mu))*(landac/(4*pi*d0))^2;
 r(k,1)=sqrt(0.5)*sqrt(variance2)*(randn(1,1)+j*randn(1,1));

 end

 for k=1:N

 variance3=((distance2(k)/d0)^(-mu))*(landac/(4*pi*d0))^2;
 h(k,1)=sqrt(0.5)*sqrt(variance3)*(randn(1,1)+j*randn(1,1));
 end

 snrr1=(abs(r).^2)*Ps/(Pn*2);

 vect_outage=snrr1>=(2^(2*R)-1);

 Nslots=Nslots+1;
 if(sum(vect_outage)>=1)

 %We look at the selection errors for 4 Relays 00,01,10,11
 % We send two pairs of bits for each 1

 vect_sec_link=abs(h).^2.*vect_outage;
 [maxim,pos1]=max(vect_sec_link); %We take the maximum
 pos1=pos1-1; % To start positions at 0
 b1=dec2bin(pos1); % We pass the positions to binary

 a=rand(1);
 b=rand(1);
 c=rand(1);
 dd=rand(1);

 l=0;

 l(4)=b1(length(b1));

 if pos1>1

 l(3)=b1(length(b1)-1);
 end

 if pos1>3
 l(2)=b1(length(b1)-2);
 end

 if pos1>7
 l(1)=b1(length(b1)-3);
 end

 l=l>'0';



 if a<=Probe
 l(4)=not(l(4));

 end

 if b<=Probe && pos1>1
 l(3)=not(l(3));
 end

 if c<=Probe && pos1>3
 l(2)=not(l(2));
 end

 if dd<=Probe && pos1>7
 l(1)=not(l(1));
 end


 %Convert binary to decimal

 choice1=l(4)+2*l(3)+4*l(2)+8*l(1);


 choice1=choice1+1;


 if N<choice1
 vect_sec_link(choice1)=0;
 end


 snrr2=(Ps/(Pn*2))*(vect_sec_link(choice1));

 channel_outage_Best=channel_outage_Best+(snrr2<(2^(2*R)-1));

 if snrr2>=(2^(2*R)-1)
 rate_mean_Best=rate_mean_Best+rateslot;
 else
 rate_mean_Best=rate_mean_Best+0;
 end

 else

 channel_outage_Best=channel_outage_Best+1;
 rate_mean_Best=rate_mean_Best+0;

 end
 end
 vector_snrm(snrm)=snrm;
 prob_outage_Best(snrm)=channel_outage_Best/I;
 prob_error_paquet_Best(snrm)=1-(1-prob_outage_Best(snrm))^(L);
 rate_relay_Best(snrm)=rate_mean_Best/Nslots;
end

for snrm=1:30
 Nslots=0;
 rate_mean_Direct=0;
 channel_outage_Direct=0;
 snrml=10^(snrm/10);
 rateslot=R;

 variance1(snrm)=((d/d0)^(-mu))*(landac/(4*pi*d0))^2;
 Ps=snrml*Pn/variance1(snrm);

 for aux=1:I

 e=sqrt(0.5)*sqrt(variance1)*(randn(1,1)+j*randn(1,1));

 snrr1=(abs(e).^2)*Ps/(Pn);

 vect_outage=snrr1>=(2^(2*R)-1);

 Nslots=Nslots+1;
 if(vect_outage==1)



 if snrr1>=(2^(2*R)-1)
 rate_mean_Direct=rate_mean_Direct+rateslot;
 else
 rate_mean_Direct=rate_mean_Direct+0;
 end

 else

 channel_outage_Direct=channel_outage_Direct+1;
 rate_mean_Direct=rate_mean_Direct+0;

 end
 end
 vector_snrm(snrm)=snrm;
 prob_outage_Direct(snrm)=channel_outage_Direct/I;
 prob_error_paquet_Direct(snrm)=1-(1-prob_outage_Direct(snrm))^(L);
 rate_relay_Direct(snrm)=rate_mean_Direct/Nslots;
end

semilogy(vector_snrm,prob_outage_Alamouti,'b-s');
hold on
grid on
%plot(vector_snrm,rate_Alamouti,'r');
semilogy(vector_snrm,prob_outage_Best,'r-*');
%plot(vector_snrm,rate_relay_Best,'k');
hold on
grid on 
semilogy(vector_snrm,prob_outage_Direct,'g-o');
hold on
%plot(vector_snrm,rate_relay_Direct,'k--');
xlabel('SNR dB');
ylabel('outage probability');
grid on
title('\fontname{Arial}outage probability without feedback error ','FontSize',12)
legend('Alamouti','Best relay','Direct')
