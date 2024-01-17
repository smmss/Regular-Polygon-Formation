close all
clear
clc
%% Initialize
% global my_phi
% simulation parameters

iterations = 1200;%120(sec)
time_step = 0.1;
tf = iterations*time_step; %final time
ts = time_step;
time = 0:ts:tf-0.1;
Tf=iterations;

Na = 8;%Number of Agents
d_L = 5;
d = sin(pi/(Na))*d_L*2;% Equation (14),  Desired Distance Between Agents
k = 1.3;
K = 1.3;
r = k*d;% Interaction range between agents
e = 0.3;% Epsilon in Equation 6
h = 0.1;% h in Equation 11
ra   = sigmanorm(r,e);
da   = sigmanorm(d,e);


%%% Parameters of Equation 4
c1a = 20;
c2a = 6;
c1g = 12;
c2g = 9;
da1_L = sigmanorm(d_L,e);
% define gamma agents

Ng = 1;% Number of leader

qg = 0*ones(Ng , 2 , iterations);%Position of leader
pg = zeros(Ng , 2 , iterations);%Velocity of leader

pg(1,:,:) = 0.5*zeros(1 , 2 , iterations);%Initialize the leader position


start_L = 1;
pg(1,:,start_L:end) = 2*ones(1 , 2 , iterations-start_L+1);
% define alpha agents


qa = zeros(Na , 2 , iterations); %Position of Agents
qa(:,:,1) = 10*(randn(Na , 2 ));
pa = zeros(Na , 2 , iterations);%Velocity of Agents
Aa = randi(Ng , Na , 1);% Assign agents to leaders

N_link = 0.5*Na*(Na-1);% Maximum number of edges between agents
distance_qa = zeros(N_link , iterations);% Maximum number of distances between agents
distance_AG = zeros(Na,iterations); % Leader & Agents
distance_AG1 = zeros(Na,iterations); % Leader & Agents

for i = 1:Na
distance_AG(i,1) = i;
end


uc = zeros(Na , 2 , iterations);
dis_A = zeros(1,iterations);
%% simulation

deltaq = zeros(Na , 2); % qa - qg
deltap = zeros(Na , 2); % pa - pg
it = 1;
num_flt = Na;
num_flt2 = 0;
time_flt = randi([150 200]);
time_flt2=time_flt;
faul = [];
faul2= [];
cnt = 1;
cnt2=0;
ts_flt = [];
tf_flt = time_flt;
my_dis1 = zeros(N_link,iterations);% Initialize all distances between agents

while it <= iterations
    %%% Fault
    if it == time_flt
        while( true )
            if num_flt <= 3
                break
            end
            rand_i = randi([1 Na]); % creat stochastic fault time
            rf = ismember(rand_i,faul);  % Checking for multiple faults not occurring for an agent
            if rf == 0
                faul = [faul rand_i];  
                cnt2 = cnt2 +1;
                break
            end
        end       
        if num_flt ~= 3
            ts_flt = [ts_flt randi([100 200])]; % step time of next faults
            num_flt = num_flt - 1;
            num_flt2 = num_flt2 + 1;
            time_flt = time_flt + ts_flt(cnt);
            tf_flt = [tf_flt time_flt]; %Times when faults occur
            cnt = cnt +1;
        end
    end
    %%%

    % calculate qdelta and pdelta
    for i = 1:(Na)
        deltaq(i,:) = qa(i,:,it) - qg(Aa(i),:,it);
        deltap(i,:) = pa(i,:,it) - pg(Aa(i),:,it);
    end
    %%%
    d = sin(pi/(Na-size(faul,2)))*d_L*2; % 
    da=sigmanorm(d,e);
    r = k*d;
    ra   = sigmanorm(r,e);
    dis_A(1,it) = d;
    %%%
    for i = 1:(Na)
         
        if any(i == faul)
            uc(i,:,it) = 0;            
        else
            uc(i,:,it) = u(i , qa(:,:,it) , pa(:,:,it) ,...
                deltaq , deltap ,h  , d_L , d, da , da1_L , ra , e ,...
                c1a , c2a  , c1g , c2g );
            pa(i,:,it+1) = pa(i,:,it) + time_step*(uc(i,:,it));
        end
    end
    %%% update
    qg(:,:,it+1) = qg(:,:,it) + time_step*pg(:,:,it);
%     pg(:,:,it+1) = pg(:,:,it);
    qa(:,:,it+1) = qa(:,:,it) + time_step*pa(:,:,it);
    %%% Total Distance Calculate
    number_link = 1;
    for i=1:(Na)
        for j=1:(Na)
            if i<j
                distance_qa(number_link,it) = norm(qa(i,:,it)-qa(j,:,it));
                number_link = number_link+1;
            end
        end
        distance_AG(i,it) = norm(qa(i,:,it)-qg(1,:,it));
    end
        my_dis1(:,it) = sortrows(distance_qa(:,it));
    %%%
    it = it + 1;
end
cnt = 1;
%%% Distance between Leader and Followers
distance_AG = sortrows(distance_AG,iterations);
%% plot Distance

%%%
tf_flt1 = [1 tf_flt iterations];
%%%
AA_Error = zeros(Na,iterations);
AL_Error = zeros(Na,iterations);
figure;
for i=1:(numel(tf_flt1)-1)    
    plot(time(:,tf_flt1(i):tf_flt1(i+1)),distance_AG(1:Na-i+1,tf_flt1(i):tf_flt1(i+1)));
    title('Distances Between Leader \& Agents','interpreter','latex','FontSize',14 )
    xlabel('Time $(s)$','interpreter','latex','FontSize',14 )
    ylabel('Distance $(m)$','interpreter','latex','FontSize',14 )
    axis ([0 tf ,0 11])
    hold on
    %%%
    AL_Error(1:Na-i+1,tf_flt1(i):tf_flt1(i+1)) =  d_L - distance_AG(1:Na-i+1,tf_flt1(i):tf_flt1(i+1));
end
print('AL_P','-depsc2','-r600')
figure;
for i=1:(numel(tf_flt1)-1)    
    plot(time(:,tf_flt1(i):tf_flt1(i+1)),my_dis1(1:Na-i+1,tf_flt1(i):tf_flt1(i+1)));
    title('Shortest Distances Between Agents','interpreter','latex','FontSize',14 )
    xlabel('Time $(s)$','interpreter','latex','FontSize',14 )
    ylabel('Distance $(m)$','interpreter','latex','FontSize',14)
    axis ([0 tf ,0 11])  
    hold on
    %%% Error Calculation
    AA_Error(1:Na-i+1,tf_flt1(i):tf_flt1(i+1)) = dis_A(:,tf_flt1(i):tf_flt1(i+1)) - my_dis1(1:Na-i+1,tf_flt1(i):tf_flt1(i+1));
end
p1 = plot(time,dis_A,'k--','DisplayName','Desired Distance');
legend(p1,'interpreter','latex')
axis ([0 tf ,0 10])
print('AA_P','-depsc2','-r600')
AL_Error = rms(AL_Error);
figure;
P2 = plot(time,AL_Error);
title('RMSE of the Distances Between Leader \& Agents','interpreter','latex','FontSize',14 )
xlabel('Time $(s)$','interpreter','latex','FontSize',14 )
ylabel('Value $(m)$','interpreter','latex','FontSize',14 )
print('AL_Error_P','-depsc2','-r600')
saveas(P2, 'AL_Error_P.png')
starter_error = iterations - 199;
figure;
P3 = plot(time(:,starter_error:end),AL_Error(:,starter_error:end));
print('AL_zoom_P','-depsc2','-r600')
saveas(P3, 'AL_zoom_P.png')
AA_Error = rms(AA_Error);
figure;
P4 = plot(time,AA_Error);
title('RMSE of the Shortest Distances Between Agents','interpreter','latex','FontSize',14 )
xlabel('Time $(s)$','interpreter','latex','FontSize',14 )
ylabel('Value $(m)$','interpreter','latex','FontSize',14 )
print('AA_Error_P','-depsc2','-r600')
saveas(P4, 'AA_Error_P.png')
figure;
P5 = plot(time(:,starter_error:end),AA_Error(:,starter_error:end));
print('AA_zoom_P','-depsc2','-r600')
saveas(P5, 'AA_zoom_P.png')
%% Final Plot
figure;
% subplot(3,1,1)
circle_x = cosd(1:360);
circle_y = sind(1:360);
Agents = scatter(qa(:,1,Tf) , qa(:,2,Tf) , 'bo');
hold on
Leaders = scatter(qg(:,1,Tf) , qg(:,2,Tf) , 'r*');
Alpha_Edge = [];
for i = 1:(Na)
    for j = i+1:(Na)
        Alpha_Edge = [Alpha_Edge plot([qa(i,1,Tf) , qa(j,1,Tf)] , [qa(i,2,Tf) , qa(j,2,Tf)])]; 
        if a_ij(qa(:,:,Tf) , i , j , e , ra , h) == 0
            Alpha_Edge(end).Color = 'none';
        else
            Alpha_Edge(end).Color = 'm';
        end
    end
end
% Gamma_Edge = [];
% for i = 1:Na
%     Gamma_Edge = [Gamma_Edge plot([qa(i,1,Tf) , qg(Aa(i),1,Tf)] , [qa(i,2,Tf) , qg(Aa(i),2,Tf)] , 'g-')]; 
% end

T = title('t = 0.00 (sec) , N = 0 , d = 0 , d_L = 0 , c^a_1 = 0 , c^a_2 = 0 , c_1^g = 0 , c_2^g = 0 , e = 0 , h = 0 ');
T.String = sprintf("t = %2.2f (sec) , N = %d , d = %d , d_L = %d , c^a_1 = %1.1f , c^a_2 = %1.1f ,  c_1^g = %2.2f , c_2^g = %2.2f , e = %1.1f , h = %1.1f " ...
                    , time_step*(Tf-1) , Na , d  , d_L , c1a , c2a , c1g , c2g , e , h );
xlabel('X','interpreter','latex','FontSize',14)
ylabel('Y','interpreter','latex','FontSize',14)
axis equal

%% Animation Plot
        d = sin(pi/(Na))*d_L*2; %
        da=sigmanorm(d,e);
        r = k*d
        ra   = sigmanorm(r,e);


f = figure;
circle_x = cosd(1:360);
circle_y = sind(1:360);

Agents = scatter(qa(:,1,1) , qa(:,2,1) , 'bo');
hold on
str_a = string(1:Na);
txt = textscatter(0.5+qa(:,1,1),0.5+qa(:,2,1),str_a);
Leaders = scatter(qg(:,1,1) , qg(:,2,1) , 'r*');

Alpha_Edge = [];
for i = 1:(Na)
    for j = i+1:(Na)
        Alpha_Edge = [Alpha_Edge plot([qa(i,1,1) , qa(j,1,1)] , [qa(i,2,1) , qa(j,2,1)])]; 
        if a_ij(qa(:,:,1) , i , j , e , ra , h) == 0
            Alpha_Edge(end).Color = 'none';
        else
            Alpha_Edge(end).Color = 'm';
        end
    end
end

% Gamma_Edge = [];
% for i = 1:(Na)
%     Gamma_Edge = [Gamma_Edge plot([qa(i,1,1) , qg(Aa(i),1,1)] , [qa(i,2,1) , qg(Aa(i),2,1)] , 'g-')]; 
% end
circles = [];
for i = 1:numel(d_L)
    circles = [circles plot(qg(1,1,1) + d_L(i)*circle_x , qg(1,2,1) + d_L(i)*circle_y , 'k--','LineWidth',0.25)];
end
%%%
T = title('$t = 0.00$ (sec)','interpreter','latex','FontSize',14);
xlabel('$x$','interpreter','latex','FontSize',14)
ylabel('$y$','interpreter','latex','FontSize',14)
axis equal

im = frame2im(getframe(f));
[X,Map] = rgb2ind(im,256);
imwrite(X , Map , 'polygon.gif' , 'gif' , 'Loopcount',inf , 'Delay' , time_step);
%%
disp(['Times when a failure occurs= ' num2str(tf_flt/10) ]);
disp(['Number of failures: ' num2str(size(faul,2)) ' , The number of agents that become faulty :' num2str(faul)]);


for i = 1:iterations
    if i ==1
%     circles = [];
%     circles2= [];
        print('pol0','-depsc2','-r600')
    end   
    T.String = sprintf("t = %2.2f (sec)" , time_step*(i));

    Agents.XData = qa(:,1,i);
    Agents.YData = qa(:,2,i);

    txt.XData = qa(:,1,i)+0.5;
    txt.YData = qa(:,2,i)+0.5;
    
    Leaders.XData = qg(:,1,i);
    Leaders.YData = qg(:,2,i);
    if i>20
    axis([qg(:,1,i)-d_L-2.5 qg(:,1,i)+d_L+2.5,qg(:,2,i)-d_L-2.5 qg(:,2,i)+d_L+2.5])
    end
%     if i>(tf_flt(1)-10)
%         axis([qg(:,1,i)-d_L-2 qg(:,1,i)+d_L+2,qg(:,2,i)-d_L-2 qg(:,2,i)+d_L+2])
%     end


    it_tmp = 1;
    for j = 1:(Na)
        for k = j+1:(Na)
            Alpha_Edge(it_tmp).XData = [qa(j,1,i) , qa(k,1,i)]; 
            Alpha_Edge(it_tmp).YData = [qa(j,2,i) , qa(k,2,i)];
            if a_ij(qa(:,:,i) , j , k , e , ra , h) == 0
                Alpha_Edge(it_tmp).Color = 'none';
            else
                Alpha_Edge(it_tmp).Color = 'm';
            end

            it_tmp = it_tmp + 1;
        end
    end
    %%% remove gamma_edge
    if i == time_flt2
        faul2 = faul(1:cnt);
        %%%
        d = sin(pi/(Na-size(faul2,2)))*d_L*2; %
        da=sigmanorm(d,e);
        r = K*d
        ra   = sigmanorm(r,e);
        %%%
        if size(faul,2) ~= cnt
            time_flt2=time_flt2+ts_flt(cnt);
            cnt = cnt+1;
        end
    end
    for j = 1:numel(d_L)
        circles(j).XData = qg(1,1,i) + d_L(j)*circle_x;
        circles(j).YData = qg(1,2,i) + d_L(j)*circle_y;       
    end  
%     for j = 1:(Na)
%         if any(j == faul2)
%             Gamma_Edge(j).XData = [];
%             Gamma_Edge(j).YData = [];
%         else
%             Gamma_Edge(j).XData = [qa(j,1,i) , qg(Aa(j),1,i)];
%             Gamma_Edge(j).YData = [qa(j,2,i) , qg(Aa(j),2,i)];
%         end
%     end
    %%
%     plot(time(i),my_dis(:,i),'k.','LineWidth',0.1);
%If you don't want to save output as GIF or EPS, you can comment the following commands       
    drawnow
    filename = strcat('pol',num2str(i));
    print(filename,'-depsc2','-r600')
    im = frame2im(getframe(f));
    [X,Map] = rgb2ind(im,256);
    imwrite(X , Map , 'polygon.gif' , 'gif' , 'WriteMode','append' , 'Delay' , time_step);
end

%% Functions
%%% Equation 6
function y = sigmanorm(z , e) 
    y = (sqrt(1 + e*norm(z)^2)-1)/e;
end
%%% Equation 10
function y = a_ij(q , i , j , e , r , h) 
%     r_alpha = sigmanorm(r,e);
    y = rho(sigmanorm(q(j,:) - q(i,:) , e)/r , h);
end
function y = sigmaepsilon(z , e) 
    y = z./sqrt(1 + e*(norm(z))^2);
end
%%% Equation 9
function y = n_ij(q , i , j , e) 
    y = sigmaepsilon(q(j,:)-q(i,:) , e);
    %y = (q(j,:)-q(i,:))/(sqrt(1+e*norm(q(j,:)-q(i,:))^2));
end
%%% Equation 9
function y = n_ri(q  , e) 
    y = sigmaepsilon(q , e);
    %y = (q /(sqrt(1+e*norm(q)^2)));
end
%%% Equation 11
function r1 = rho(z , h)
    r1 = zeros(size(z));
    for i = 1:numel(z)
        if z(i) < h
            r1(i) = 1;
        elseif z(i) < 1
            r1(i) = 0.5*(1 + cos(pi*(z(i) - h)/(1 - h)));
        end
    end
end
%%% Equation 7
function y = phi_alpha(z , h , d, da , r)
    aa = 5;
    bb = 5;
    cc = abs(aa - bb)/sqrt(4*aa*bb);
    phi = @(z) 0.5*((aa+bb)*sigma_11(z/d+cc)+(aa-bb));
    y = rho(z/r , h).*phi(z - da);
end
%%% Equation 7
function y = phi_alpha_L(z  ,d, d1)
    aa = 5;
    bb = 5;
    cc = abs(aa - bb)/sqrt(4*aa*bb);
    phi = @(z) 0.5*((aa+bb)*sigma_11(z/d+cc)+(aa-bb));
        y = 1*phi(z - d1);
end
%%% Equation 8
function y = sigma_11(z)
    y = z/(1 + (z)^2)^0.5;
end
%%% Equation 4
function y = u(i , q , p , deltaq , deltap , h  , d_L , d, da , da1_L , ra  , e  , c1a , c2a  , c1g , c2g )
    y = [0,0];
    Na = size(q , 1);
    for j = 1:Na
        if j ~= i
            y = y + c1a*phi_alpha(sigmanorm(q(j,:)-q(i,:),e) , h , d, da , ra).*n_ij(q , i , j , e)...
                  + c2a*a_ij(q , i , j , e , ra , h)*(p(j,:)-p(i,:)); 
        end
    end
    y = y - c1g*phi_alpha_L(sigmanorm(deltaq(i,:),e)  ,d_L , da1_L )... 
        .*n_ri(deltaq(i,:) , e) - c2g*deltap(i,:);
end
