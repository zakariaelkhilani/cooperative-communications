%close all;
clear all;
I=30000;
N=10;
R=2;
d=100;
par=d/(2*(N-1));
if rem(N,2)==0
 limite=N/2;
end
if rem(N,2)==1
 limite=N/2+0.5;
end
for aux=1:limite
 distance(aux)=sqrt((d/4+par*(aux-1))^2+(d/4)^2);
end
for aux=(limite+1):N
 distance(aux)=sqrt((d/4+par*(aux-limite-1))^2+(d/4)^2);
end
distance(N)=sqrt((d/4+d/2)^2+(d/4)^2);
distance(limite)=sqrt((d/4+d/2)^2+(d/4)^2);
for aux=1:limite
 distance2(aux)=distance(limite+1-aux);
end
for aux=(limite+1):N
 distance2(aux)=distance(N+1+limite-aux);
end
fc=2.5*10^9;
c=3*10^8;
landac=c/fc;
d0=1;
mu=3;
Pn=1;
for snrm=1:30
 channel_outage=0;
 snrml=10^(snrm/10);

 variance1(snrm)=((d/d0)^(-mu))*(landac/(4*pi*d0))^2;
 Ps=snrml*Pn/variance1(snrm);
 for aux=1:I
 for cont=1:N

 variance2=((distance(cont)/d0)^(-mu))*(landac/(4*pi*d0))^2;
 r(cont,1)=sqrt(0.5)*sqrt(variance2)*(randn(1,1)+j*randn(1,1));

 end

 for cont=1:N

 variance3=((distance2(cont)/d0)^(-mu))*(landac/(4*pi*d0))^2;
 h(cont,1)=sqrt(0.5)*sqrt(variance3)*(randn(1,1)+j*randn(1,1));
 end

 snrr1=(abs(r).^2)*Ps/(Pn*2);
 vect_outage=snrr1>=(2^(2*R)-1);
 if(sum(vect_outage)>=2)

 vect_sec_link=sort(abs(h).^2.*vect_outage,'descend');
 snrr2=(Ps/(Pn*4))*(sum(vect_sec_link(1:2)));
 channel_outage=channel_outage+(snrr2<(2^(2*R)-1));

 else

 channel_outage=channel_outage+1;

 end
 end
 vector_snrm(snrm)=snrm;
 prob_outage_Alamouti(snrm)=channel_outage/I;
end
%position1=find(prob_outage<1e-3);
%valor1=min(position1);
%Ps_Alamouti=Pn*vector_snrm(valor1)/variance1(valor1);
%Direct link
for snrm=1:30
 channel_outage=0;
 snrml=10^(snrm/10);
 variance1(snrm)=((d/d0)^(-mu))*(landac/(4*pi*d0))^2;
 Ps=snrml*Pn/variance1(snrm);
for aux=1:I
r=sqrt(0.5)*sqrt(variance1(snrm))*(randn(1,1)+j*randn(1,1));
snrr1=(abs(r).^2)*(Ps/Pn);
if snrr1<((2^R)-1)
 channel_outage=channel_outage+1;
end
end
prob_outage_direct(snrm)=channel_outage/I;
end
 hold on
 grid on
semilogy(vector_snrm,prob_outage_Alamouti,'k->');
xlabel('SNR dB');
ylabel('Outage proability ');
axis([0 30 10^-4 1])
%legend('alamouti_rectangular')

%position2=find(prob_outage_direct<1e-3);
%valeur2=min(position2);
%Ps_direct=Pn*vector_snrm(valeur2)/variance1(valeur2);