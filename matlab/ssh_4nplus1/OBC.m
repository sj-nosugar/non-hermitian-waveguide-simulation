clear;clc;tic;
% SSH model obc
%% 
N = 80;
t = 1;          %hopping t 
de = -0.5   ;    %hopping delta
 t + de
  t - de
%% ham

H = zeros(2*N,2*N);

for i = 1:2:2*N-2
    
    H(i,i+1) = t + de;
    H(i+1,i) = t + de;
    H(i+1,i+2) = t - de;
    H(i+2,i+1) = t - de;
    
end

H(2*N,2*N-1) = t + de;
H(2*N-1,2*N) = t + de;

[v,e] = eig(H);

%% plot
figure(1)
plot(1:2*N,diag(e),'k.','MarkerSize',8)
%grid on
xlim([1 2*N])

figure
pM = [1:20,2*N-19:2*N];
plot(1:160,-v(:,N+40),'r*-','MarkerSize',4,'LineWidth',1)
ylim([-0.6 0.6])
e(N+40,N+40)
figure
pM = [1:20,2*N-19:2*N];
plot(1:40,-v(pM,N),'r*-','MarkerSize',4,'LineWidth',1)
hold on
plot(1:40,-v(pM,N+1),'bo-','MarkerSize',4,'LineWidth',1)
view([0,90])
%hc=log(JR)-log(V);
grid on
hold off
ylim([-0.6 1])
xlim([1 40+0.1])
set(gcf,'unit','centimeters','position',[10 5 9 7])
xlabel('$h$','Interpreter',"latex",'FontSize',14,'FontWeight','bold','Color','k')
ylabel('$|E|$','Interpreter',"latex",'FontSize',14,'FontWeight','normal','Color','k')
zlabel('1/(IPR)','Interpreter',"tex",'FontSize',10,'FontWeight','normal','Color','k')
text(0,0,'(a)','Interpreter',"none",'VerticalAlignment',"middle",'FontSize',14,"HorizontalAlignment","left")
grid off
colorbar



%%
toc;