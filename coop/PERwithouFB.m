clear all;
%close all;
I=300000;
N=7;
%R=2;
d=100;
%Calculate the distncies of each destination.
distances=d*rand(N,1);
distancey=(d/2)*rand(N,1);
distance=sqrt(distances.^2+distancey.^2);
distancies2=d-distances;
distance2=sqrt(distancies2.^2+distancey.^2);
fc=2.5*10^9;
c=3*10^8;
landac=c/fc;
d0=1;
mu=3;
Pn=1;
L=10; %Number of symbols in a package.
%PER Direct
for snrm=1:30

 rate_mean=0;
 fed_error=0;
 Nslots=0;
 error_symbol_Direct=0;
 error_symbol_Best=0;
 error_symbol_Alamouti=0;
 snrml=10^(snrm/10);
 Pesym=10^(-1);

 variance1=((d/d0)^(-mu))*(landac/(4*pi*d0))^2;
 Ps=snrml*Pn/variance1;

 for aux=1:I

 e=sqrt(0.5)*sqrt(variance1)*(randn(1,1)+j*randn(1,1));

 snr_direct=(abs(e).^2)*Ps/(Pn);
 Pesym_direct=2*(1/2*erfc(sqrt(snr_direct/2)))-(1/2*erfc(sqrt(snr_direct/2))).^2;
 outage_direct=Pesym_direct<Pesym;

 if outage_direct==0
 error_symbol_Direct=error_symbol_Direct+1;
 else
 error_symbol_Direct=error_symbol_Direct+0;
 end
 end

 vector_snrm(snrm)=snrm;
 prob_outage_Direct(snrm)=error_symbol_Direct/I;
 prob_error_paquet_Direct(snrm)=1-(1-prob_outage_Direct(snrm))^(L);
 %%Best Relay



 for aux=1:I

 for k=1:N
 variance2=((distance(k)/d0)^(-mu))*(landac/(4*pi*d0))^2;
 r(k,1)=sqrt(0.5)*sqrt(variance2)*(randn(1,1)+j*randn(1,1));
 end

 for k=1:N
 variance3=((distance2(k)/d0)^(-mu))*(landac/(4*pi*d0))^2;
 h(k,1)=sqrt(0.5)*sqrt(variance3)*(randn(1,1)+j*randn(1,1));
 end


 %We calculated the decoding set

 snr_relays=(abs(r).^2)*Ps/(Pn*2);
 Pesym_relays=2*(1/2*erfc(sqrt(snr_relays/2)))-(1/2*erfc(sqrt(snr_relays/2))).^2;
 outage_relays=Pesym_relays<Pesym;
 Pesym_relays=Pesym_relays.*outage_relays;

 if sum(outage_relays(1:end))==0
 Pesym_relays=inf;
 else
 Pesym_relays=sum(Pesym_relays(1:end))/sum(outage_relays(1:end));
 end



 %We look at the relays that are in outage to transmit to the Destination
 vect_sec_link=abs(h).^2.*outage_relays;
 snr_relays2=(Ps/(Pn*2)).*vect_sec_link;
 Pesym_relays2=2*(1/2*erfc(sqrt(snr_relays2/2)))-(1/2*erfc(sqrt(snr_relays2/2))).^2;
 Pesym_relays2=Pesym_relays2.*outage_relays;
 outage_relays2=Pesym_relays2<Pesym;

 if sum(outage_relays2(1:end))==0
 Pesym_relays2=inf;
 else
 Pesym_relays2=sum(Pesym_relays2(1:end))/sum(outage_relays2(1:end));
 end


 Pebit=Pesym_relays2/2;

 %We calculate the best positions

 vect_sec_link=abs(h).^2.*outage_relays;
 [maxim,pos1]=max(vect_sec_link); %We take the maximum
 vect_sec_link2=vect_sec_link(pos1); %We keep it
 vect_sec_link(pos1)=0; %Put your position to 0
 [maxim2,pos2]=max(vect_sec_link); %We take the next maximum
 vect_sec_link(pos1)=vect_sec_link2;

 a=rand(1);
 b=rand(1);
 c=rand(1);
 dd=rand(1);


 pos1=pos1-1; %To start positions at 0
 pos2=pos2-1;
 b1=dec2bin(pos1); %We pass the positions to binary
 b2=dec2bin(pos2);
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
 if a<=Pebit
 l(4)=not(l(4));
 end
 if b<=Pebit && pos1>1
 l(3)=not(l(3));
 end
 if c<=Pebit && pos1>3
 l(2)=not(l(2));
 end
 if dd<=Pebit && pos1>7
 l(1)=not(l(1));
 end
 %Convert binary to decimal
 choice1=l(4)+2*l(3)+4*l(2)+8*l(1);
 choice1=choice1+1;

 if N<choice1
 vect_sec_link(choice1)=0;
 end


 if(sum(outage_relays)>0)


 snr_Best=(Ps/(Pn*2))*vect_sec_link(choice1);
 Pesym_Best=2*(1/2*erfc(sqrt(snr_Best/2)))-(1/2*erfc(sqrt(snr_Best/2))).^2;
 outage_Best=Pesym_Best<Pesym;
 if outage_Best==0
 error_symbol_Best=error_symbol_Best+1;
 else
 error_symbol_Best=error_symbol_Best+0;
 end
 else
 error_symbol_Best=error_symbol_Best+1;
 rate_mean=rate_mean+0;
 end
 end

 
 vector_snrm(snrm)=snrm;
 prob_outage_Best(snrm)=error_symbol_Best/I;
 prob_error_paquet_Best(snrm)=1-(1-prob_outage_Best(snrm))^(L);



 %Alamouti

 for aux=1:I

 for k=1:N
 variance2=((distance(k)/d0)^(-mu))*(landac/(4*pi*d0))^2;
 r(k,1)=sqrt(0.5)*sqrt(variance2)*(randn(1,1)+j*randn(1,1));
 end

 for k=1:N
 variance3=((distance2(k)/d0)^(-mu))*(landac/(4*pi*d0))^2;
 h(k,1)=sqrt(0.5)*sqrt(variance3)*(randn(1,1)+j*randn(1,1));
 end
 %We calculated the decoding set

 snr_relays=(abs(r).^2)*Ps/(Pn*2);
 Pesym_relays=2*(1/2*erfc(sqrt(snr_relays/2)))-(1/2*erfc(sqrt(snr_relays/2))).^2;
 outage_relays=Pesym_relays<Pesym;
 Pesym_relays=Pesym_relays.*outage_relays;

 if sum(outage_relays(1:end))==0
 Pesym_relays=inf;
 else
 Pesym_relays=sum(Pesym_relays(1:end))/sum(outage_relays(1:end));
 end

 %We look at the relays that are in outage to transmit to the Destination

 vect_sec_link=abs(h).^2.*outage_relays;
 snr_relays2=(Ps/(Pn*2)).*vect_sec_link;
 Pesym_relays2=2*(1/2*erfc(sqrt(snr_relays2/2)))-(1/2*erfc(sqrt(snr_relays2/2))).^2;
 Pesym_relays2=Pesym_relays2.*outage_relays;
 outage_relays2=Pesym_relays2<Pesym;


 if sum(outage_relays2(1:end))==0
 Pesym_relays2=inf;
 else
 Pesym_relays2=sum(Pesym_relays2(1:end))/sum(outage_relays2(1:end));
 end

 Pebit=Pesym_relays2/2;

 %We calculate the best positions

 vect_sec_link=abs(h).^2.*outage_relays;
 [maxim,pos1]=max(vect_sec_link); %take the maximum
 vect_sec_link2=vect_sec_link(pos1); %We keep it
 vect_sec_link(pos1)=0; %position a 0
 [maxim2,pos2]=max(vect_sec_link); %second maximum
 vect_sec_link(pos1)=vect_sec_link2;

 a=rand(1);
 b=rand(1);
 c=rand(1);
 dd=rand(1);

 pos1=pos1-1; %To start positions at 0
 pos2=pos2-1;
 b1=dec2bin(pos1); %We pass the positions to binary

 b2=dec2bin(pos2);
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
 if a<=Pebit
 l(4)=not(l(4));
 end
 if b<=Pebit && pos1>1
 l(3)=not(l(3));
 end
 if c<=Pebit && pos1>3
 l(2)=not(l(2));
 end
 if dd<=Pebit && pos1>7
 l(1)=not(l(1));
 end
 %convert binary to decimal
 choice1=l(4)+2*l(3)+4*l(2)+8*l(1);
 %Now we do the same for the second selection
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
 a=rand(1);
 b=rand(1);
 c=rand(1);
 dd=rand(1);
 if a<=Pebit
 l(4)=not(l(4));
 end
 if b<=Pebit && pos1>1
 l(3)=not(l(3));
 end
 if c<=Pebit && pos1>3
 l(2)=not(l(2));
 end
 if dd<=Pebit && pos1>7
 l(1)=not(l(1));
 end
 %Convertion binary to decimal
 choice2=l(4)+2*l(3)+4*l(2)+8*l(1);
 choice1=choice1+1;
 choice2=choice2+1;

 if N<choice1
 vect_sec_link(choice1)=0;
 end
 if N<choice2
 vect_sec_link(choice2)=0;
 end

 if(sum(outage_relays)>=2)

 snr_Alamouti=(Ps/(Pn*4))*(vect_sec_link(choice1)+vect_sec_link(choice2));
 Pesym_Alamouti=2*(2*(1/2*erfc(sqrt(snr_Alamouti/2)))-(1/2*erfc(sqrt(snr_Alamouti/2))).^2);
 outage_Alamouti=Pesym_Alamouti<Pesym;
 if outage_Alamouti==0
 error_symbol_Alamouti=error_symbol_Alamouti+1;
 else
 error_symbol_Alamouti=error_symbol_Alamouti+0;
 end
 else
 error_symbol_Alamouti=error_symbol_Alamouti+1;
 end
 end

 vector_snrm(snrm)=snrm;
 prob_outage_Alamouti(snrm)=error_symbol_Alamouti/I;
 prob_error_paquet_Alamouti(snrm)=1-(1-prob_outage_Alamouti(snrm))^(L);
end
semilogy(vector_snrm,prob_error_paquet_Direct,'g-o');
hold on
grid on
%plot(vector_snrm,rate_Alamouti,'r');
semilogy(vector_snrm,prob_error_paquet_Best,'r-*');
%plot(vector_snrm,rate_relay_Best,'k');
hold on
grid on 
semilogy(vector_snrm,prob_error_paquet_Alamouti,'b-s');
hold on
%plot(vector_snrm,rate_relay_Direct,'k--');
xlabel('SNR dB');
ylabel('PER');
grid on
%title('\fontname{Arial}PER using 10relays without feedback error ','FontSize',12)
legend('Direct','Best relay','Alamouti')