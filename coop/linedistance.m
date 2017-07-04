clear all;
I=30000;
N=10;
R=2;
d=100;
par=d/(2*(N-1));
 for aux=1:N
 distances(aux)=d/4+par*(aux-1);
 end
 for aux=1:N
 distances(aux)=d/4+par*(aux-1);
 end


for aux=1:N
 distance2(aux)=distances(N+1-aux);
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

 variance2=((distances(cont)/d0)^(-mu))*(landac/(4*pi*d0))^2;
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
 prob_outage_Alamoutiline(snrm)=channel_outage/I;
end
%posicion1=find(prob_outage_Alamoutiline<1e-3);
%valor1=min(posicion1);
%Ps_Alamouti=Pn*vector_snrm(valor1)/variance1(valor1);
%Diract link
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
semilogy(vector_snrm,prob_outage_direct,'b-o');
hold on
grid on
semilogy(vector_snrm,prob_outage_Alamoutiline,'k-o');
 xlabel('SNR dB')
 ylabel('Outage Probability');
 axis([0 30 10^-4 1])
 %legend('Direct link','Alamouti line')
%posicion2=find(prob_outage_direct<1e-3);
%valor2=min(posicion2);
%Ps_directo=Pn*vector_snrm(valor2)/variance1(valor2);