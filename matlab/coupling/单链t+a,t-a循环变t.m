clear all

a=3/4;
v=1;
Nc2=40;
tset=(-3:0.01:3);
ra0Set = zeros(Nc2, numel(tset));
ra1Set = zeros(Nc2, numel(tset));
ra2Set = zeros(Nc2, numel(tset));
% ra3Set = zeros(Nc2, numel(t0set));


for ii = 1 : 1 : numel(tset)
t=tset(ii);  
k1=t-a;
k2=t+a;
Ham1=zeros(Nc2);
for nc1=(1:Nc2/2) 
    %(2,1)
        Ham1(2*nc1,2*nc1-1)=k1;
         %(1,2)
        Ham1(2*nc1-1,2*nc1)=k2;
end
 for nc1=(1:(Nc2-1)/2)
     %(2,3)
        Ham1(2*nc1,2*nc1+1)=v;
          %(3,2)
        Ham1(2*nc1+1,2*nc1)=v;
 end
  % Ham1(1,Nc2)=v;
% Ham1(Nc2,1)=v;

[B2,c2]=eig(Ham1);
x1=diag(c2);
ra0Set(:,ii)=x1;
x2=real(x1);
x3=imag(x1);
[x22,I]=sortrows(x2);
D(:,1)=x2;
D(:,2)=x3;
A=D(I,:);
x4=(x2.^2+x3.^2).^(0.5);
x44=sort(x4);
ra0Set(:,ii)=(D(:,1));
ra1Set(:,ii)=(D(:,2));
ra2Set(:,ii)=(x44);
end
t1matx=zeros(Nc2,numel(tset));
for ii = 1 : 1 : numel(tset)
    
     t1=tset(ii);
   
t1matx(:,ii)=t1;
end
figure
%subplot(1,2,2) % Energy bands
plot(t1matx,(ra2Set),'k.','MarkerFaceColor','[0.8 0.8 0.8]')

%plot(t1matx,(ra2Set),'.','color','[0.8 0.8 0.8] ','MarkerFaceColor','[0.8 0.8 0.8]')
% axis([-3 3 0 4]);
xlabel('t');
ylabel('|E|');


 %title(['k_1=',num2str(k1),',k_2=',num2str(k2),',k_{11}=',num2str(k11),',k_{22}=',num2str(k22),',v=',num2str(v),',N_{site}=',num2str(Nc2)]);
 %title(['k_1=',num2str(k1),',k_2=',num2str(k2),',k_3=',num2str(k22),',v=',num2str(v),',N_{site}=',num2str(Nc2)]);
 %title(['k_1=k_2=',num2str(k2),',k_3=',num2str(k22),',v=',num2str(v),',N_{site}=',num2str(Nc2)]);
 
%subplot(1,2,1) % Energy bands
figure
plot(t1matx,ra0Set,'r .','MarkerFaceColor','r')
% axis([-3 3 0 4]);
xlabel('t');
ylabel('Re|E|');
% ylim([-0.1,0.1]);
%title(['k_1=',num2str(k1),',k_2=',num2str(k2),',k_3=',num2str(k22),',v=',num2str(v),',N_{site}=',num2str(Nc2)]);
% title(['\kappa_1=',num2str(k1),'\nu',',\kappa_2=',num2str(k2),'\nu',',\kappa_3=',num2str(k22),'\nu',',N_{site}=',num2str(Nc2)]);
 title(['a=',num2str(a),',\nu=',num2str(v),',N_{site}=',num2str(Nc2)]);
