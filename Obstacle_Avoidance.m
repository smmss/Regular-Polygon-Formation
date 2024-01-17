close all
clear
clc
%% Initialize

% simulation parameters

iterations = 1200;%120(sec)
time_step = 0.1;
tf = iterations*time_step; %final time
ts = time_step;
time = 0:ts:tf-0.1;
Tf=iterations;
% constants
k = 1.3; 
e = 0.3; % Epsilon in Equation 5
h = 0.1; % h in Equation 10
d_L = 5; %Radius of circle formation
d_gamma = sigmanorm(d_L,e);

%%% Parameters of Equation 4
c1a = 10;
c2a = 4;
c1b = 40;
c2b = 10;
c1g = 7;
c2g = 5;
%%%

%%% define gamma agents

Ng = 1; % Number of leader

qg = 0*ones(Ng , 2 , iterations); %Position of leader
pg = zeros(Ng , 2 , iterations);  %Velocity of leader
start_L = 400;                      %Start time of leader to move
pg(1,:,:) = zeros(1 , 2 , iterations);%Initialize the leader position
pg(1,:,start_L:end) = 2*ones(1 , 2 , iterations-start_L+1);%Initialize the leader velocity 

%%% define alpha agents

Na = 12; %Number of Agents
qa = zeros(Na , 2 , iterations); %Position of Agents

qa(:,:,1) = 10*randn(Na , 2 );
pa = zeros(Na , 2 , iterations);%Velocity of Agents
Aa = randi(Ng , Na , 1);% Assign agents to leaders

N_link = 0.5*Na*(Na-1);% Maximum number of edges between agents
distance_qa = zeros(N_link , iterations);% Maximum number of distances between agents
distance_AG = zeros(Na,iterations); % ... Leader & Agents

% define beta agents

Nb = 3;% Number of Obstacles
yk = zeros(Nb , 2 , iterations);% For center of obstacles
Rk = [5,2.5,2.5]; %Radius of obstacles 
yk(1,:,:) = repmat([15 15], 1, 1, size(yk, 3)); %Initialize the center of obstacles
yk(2,:,:) = repmat([44 55], 1, 1, size(yk, 3)); %...center
yk(3,:,:) = repmat([55 44], 1, 1, size(yk, 3)); %...center
q_hat = zeros(Na*Nb , 2 , iterations);
p_hat = zeros(Na*Nb , 2 , iterations);
d_obs = 3; % The distance between agents and the surface of obstacles (d')
d_beta = sigmanorm(d_obs,e);


uc = zeros(Na , 2 , iterations); %Initialize the control input
%% simulation

deltaq = zeros(Na , 2); % qa - qg
deltap = zeros(Na , 2); % pa - pg
it = 1;

%%%
d = sin(pi/Na)*d_L*2; % Equation (14),  Desired Distance Between Agents
d_alpha = sigmanorm(d,e);
r = k*d; % Interaction range between agents
r_alpha= sigmanorm(r,e);
%%%
my_dis1 = zeros(N_link,iterations); % Initialize all distances between agents
dis_AG = zeros(Na,iterations); % Initialize all distances between agents & leader
while it <= iterations
    % calculate qdelta and pdelta
    for i = 1:Na
        deltaq(i,:) = qa(i,:,it) - qg(Aa(i),:,it);
        deltap(i,:) = pa(i,:,it) - pg(Aa(i),:,it);
    end

    for i = 1:Na
        for j = 1:Nb
            %%% Equation 5
            mu = Rk(j)/(norm(qa(i,:,it)-yk(j,:,it)));
            ak = (qa(i,:,it)-yk(j,:,it))/(norm(qa(i,:,it)-yk(j,:,it)));
            P = eye - ak*transpose(ak);
            q_hat((Na*j-Na+i),:,it) = mu*qa(i,:,it)+(1-mu)*yk(j,:,it);
            p_hat((Na*j-Na+i),:,it) = mu*P*pa(i,:,it);
            if norm(q_hat((Na*j-Na+i),:,it)-qa(i,:,it))<0.9
                norm(q_hat((Na*j-Na+i),:,it)-qa(i,:,it))
            end
            %%%
        end
        uc(i,:,it) = u(i , qa(:,:,it) , pa(:,:,it) ,...
                                       deltaq , deltap ,...
                                       h  , d_L , d , d_alpha , d_gamma ...
                                       , r_alpha   , e  , c1a , c2a , c1b , c2b , c1g , ...
                                       c2g, yk(:,:,it), Rk, d_beta, ...
                                       q_hat(:,:,it),p_hat(:,:,it));           
    end
    % update   
    qg(:,:,it+1) = qg(:,:,it) + time_step*pg(:,:,it);
    pa(:,:,it+1) = pa(:,:,it) + time_step*(uc(:,:,it));
    qa(:,:,it+1) = qa(:,:,it) + time_step*pa(:,:,it);
    %%% Calculate all distances
    number_link = 1;
    for i=1:Na
        for j=1:Na
            if i<j
                distance_qa(number_link,it) = norm(qa(i,:,it)-qa(j,:,it));
                number_link = number_link+1;
            end
        end
        distance_AG(i,it) = norm(qa(i,:,it)-qg(1,:,it));
    end
    my_dis1(:,it) = sortrows(distance_qa(:,it));
    dis_AG(:,it) = sortrows(distance_AG(:,it));
    %%%
    it = it + 1;
end
distance_qa(:,end) = distance_qa(:,end-1);  % Shift Data
my_dis1(:,end) = my_dis1(:,end-1);
dis_AG(:,end) = dis_AG(:,end-1);
%% plot Distance
figure;
plot(time,distance_qa);
title('All Distance Between Each Agent')
xlabel('Time')
ylabel('Distance')
% axis ([0 tf ,0 50])
% print('LA_Dis_obs','-depsc2','-r600')

figure;
plot(time,dis_AG);
title('All Distances Between Leader \& Agents','interpreter','latex','FontSize',14)
xlabel('Time $(s)$','interpreter','latex','FontSize',14)
ylabel('Distance  $(m)$','interpreter','latex','FontSize',14)
axis ([0 tf ,0 20])
print('AL_O','-depsc2','-r600')
figure;    
plot(time,my_dis1(1:Na,:));
title('Shortest Distances Between Agents','interpreter','latex','FontSize',14)
xlabel('Time $(s)$','interpreter','latex','FontSize',14)
ylabel('Distance $(m)$','interpreter','latex','FontSize',14)
axis ([0 tf ,0 9])
print('AA_O','-depsc2','-r600')
%%% Error: Agent & Leader 
Error_AL = rms(dis_AG - d_L);
figure;
P1 = plot(time,Error_AL);
title('RMSE of the Distances Between Leader \& Agents','interpreter','latex','FontSize',14)
xlabel('Time $(s)$','interpreter','latex','FontSize',14)
ylabel('Value $(m)$','interpreter','latex','FontSize',14)
axis ([0 tf ,0 15])
print('AL_Error_O','-depsc2','-r600')
saveas(P1, 'AL_Error_O.png')
%%% Error: Agent-Agent
Error_AA = rms(my_dis1(1:Na,:) - d);
figure;
P2 = plot(time,Error_AA);
title('RMSE of the Shortest Distances Between Agents','interpreter','latex','FontSize',14)
xlabel('Time $(s)$','interpreter','latex','FontSize',14)
ylabel('Value $(m)$','interpreter','latex','FontSize',14)
axis ([0 tf ,0 10])
print('AA_Error_O','-depsc2','-r600')
saveas(P2, 'AA_Error_O.png')
%%% 
steady_error = iterations - 199;
figure;
P3 = plot(time(:,steady_error:end),Error_AL(:,steady_error:end));
% axis ([0 tf ,0 10])
print('AL_Error_Zoom_O','-depsc2','-r600')
saveas(P3, 'AL_Error_Zoom_O.png')
%%% Error: Agent-Agent
figure;
P4 = plot(time(:,steady_error:end),Error_AA(:,steady_error:end));
% axis ([0 tf ,0 10])
print('AA_Error_Zoom_O','-depsc2','-r600')
saveas(P4, 'AA_Error_Zoom_O.png')
%% Final Plot
figure;
% subplot(3,1,1)
circle_x = cosd(1:360);
circle_y = sind(1:360);

Agents = scatter(qa(:,1,Tf) , qa(:,2,Tf) , 'bo');
hold on
Leaders = scatter(qg(:,1,Tf) , qg(:,2,Tf) , 'r*');
Obstacles = [];
for i = 1:Nb
    Obstacles = [Obstacles plot(yk(i,1,Tf) + Rk(i)*circle_x , yk(i,2,Tf) + Rk(i)*circle_y , 'k-')];
end
%%% Edge between agents
Alpha_Edge = [];
for i = 1:Na
    for j = i+1:Na
        Alpha_Edge = [Alpha_Edge plot([qa(i,1,Tf) , qa(j,1,Tf)] , [qa(i,2,Tf) , qa(j,2,Tf)])]; 
        if a_ij(qa(:,:,Tf) , i , j , e , r_alpha , h) == 0
            Alpha_Edge(end).Color = 'none';
        else
            Alpha_Edge(end).Color = 'm';
        end
    end
end
%%% Edge between agents and leader
% Gamma_Edge = [];
% for i = 1:Na
%     Gamma_Edge = [Gamma_Edge plot([qa(i,1,Tf) , qg(Aa(i),1,Tf)] , [qa(i,2,Tf) , qg(Aa(i),2,Tf)] , 'g-')]; 
% end

T = title('t = 0.00 (sec) , N = 0 , d = 0 , d_L = 0 , c^a_1 = 0 , c^a_2 = 0 , c_1^g = 0 , c_2^g = 0 , e = 0 , h = 0 ');
T.String = sprintf("t = %2.2f (sec) , N = %d , d = %d , d_L = %d , c^a_1 = %1.1f , c^a_2 = %1.1f ,  c_1^g = %2.2f , c_2^g = %2.2f , e = %1.1f , h = %1.1f " ...
                    , time_step*(Tf-1) , Na , d  , d_L , c1a , c2a , c1g , c2g , e , h );
xlabel('X')
ylabel('Y')
axis equal
% grid on 
% grid minor

%% Animation Plot

f = figure;
% subplot(2,1,1)

circle_x = cosd(1:360);
circle_y = sind(1:360);

Agents = scatter(qa(:,1,1) , qa(:,2,1) , 'bo');
hold on
% Number 
% str_a = string(1:Na);
% txt = textscatter(2+qa(:,1,1),2+qa(:,2,1),str_a);
Leaders = scatter(qg(:,1,1) , qg(:,2,1) , 'r*');
Obstacles = [];
for i = 1:Nb
    Obstacles = [Obstacles plot(yk(i,1,1) + Rk(i)*circle_x , yk(i,2,1) + Rk(i)*circle_y , 'k-')];
end
%%% Edge between agents
Alpha_Edge = [];
for i = 1:Na
    for j = i+1:Na
        Alpha_Edge = [Alpha_Edge plot([qa(i,1,1) , qa(j,1,1)] , [qa(i,2,1) , qa(j,2,1)])]; 
        if a_ij(qa(:,:,1) , i , j , e , r , h) == 0
            Alpha_Edge(end).Color = 'none';
        else
            Alpha_Edge(end).Color = 'm';
        end
    end
end
%%% Edge between agents and obstacles
Beta_Edge = [];
for i = 1:Na
    for j = 1:Nb
        Beta_Edge = [Beta_Edge plot([qa(i,1,1) , q_hat((Na*j-Na+i),1,1)] , [qa(i,2,1) , q_hat((Na*j-Na+i),2,1)])]; 
        if b_ik([qa(:,:,1);q_hat(:,:,1)], i , (Na*j+i) , e , d_beta , h) == 0
            Beta_Edge(end).Color = 'none';
        else
            Beta_Edge(end).Color = 'r';
        end
    end
end
Gamma_Edge = [];
%%% Edge between agents and leader
% for i = 1:Na
%     Gamma_Edge = [Gamma_Edge plot([qa(i,1,1) , qg(Aa(i),1,1)] , [qa(i,2,1) , qg(Aa(i),2,1)] , 'g-')]; 
% end
circles = [];
for i = 1:numel(d_L)
    circles = [circles plot(qg(1,1,1) + d_L(i)*circle_x , qg(1,2,1) + d_L(i)*circle_y , 'k--','LineWidth',0.25)];
end

T = title('$t = 0.00$ (sec)','interpreter','latex','FontSize',14);
xlabel('$x$','interpreter','latex','FontSize',14)
ylabel('$y$','interpreter','latex','FontSize',14)
axis equal

%%% gif
im = frame2im(getframe(f));
[X,Map] = rgb2ind(im,256);
imwrite(X , Map , 'obstacle.gif' , 'gif' , 'Loopcount',inf , 'Delay' , time_step);
%%
for i = 1:iterations
    if i ==1
%     circles = [];
%     circles2= [];
        print('obs_0','-depsc2','-r600')
    end
    if i>(start_L) && (i<start_L+200)
        axis([qg(:,1,i)-3*d_L qg(:,1,i)+3*d_L,qg(:,2,i)-3*d_L qg(:,2,i)+3*d_L])
    elseif (i>start_L+200)
        axis([qg(:,1,i)-2*d_L qg(:,1,i)+2*d_L,qg(:,2,i)-2*d_L qg(:,2,i)+2*d_L])
    end
    T.String = sprintf("t = %2.2f (sec)" , time_step*(i));

    Agents.XData = qa(:,1,i);
    Agents.YData = qa(:,2,i);

    Leaders.XData = qg(:,1,i);
    Leaders.YData = qg(:,2,i);

    for j = 1:Nb
        Obstacles(j).XData = yk(j,1,i) + Rk(j)*circle_x;
        Obstacles(j).YData = yk(j,2,i) + Rk(j)*circle_y;
    end

    it_tmp = 1;
    it_tmp1 = 1;
    for j = 1:Na
        for k = j+1:Na
            Alpha_Edge(it_tmp).XData = [qa(j,1,i) , qa(k,1,i)]; 
            Alpha_Edge(it_tmp).YData = [qa(j,2,i) , qa(k,2,i)];
            if a_ij(qa(:,:,i) , j , k , e , r_alpha , h) == 0
                Alpha_Edge(it_tmp).Color = 'none';
            else
                Alpha_Edge(it_tmp).Color = 'm';
            end

            it_tmp = it_tmp + 1;
        end
    end
    for j = 1:Na
        for k = 1:Nb
            Beta_Edge(it_tmp1).XData = [qa(j,1,i) , q_hat((Na*k-Na+j),1,i)]; 
            Beta_Edge(it_tmp1).YData = [qa(j,2,i) , q_hat((Na*k-Na+j),2,i)];
            if b_ik([qa(:,:,i);q_hat(:,:,i)] , j , (Na*k+j) , e , d_beta , h) == 0
                Beta_Edge(it_tmp1).Color = 'none';
            else
                Beta_Edge(it_tmp1).Color = 'r';
            end

            it_tmp1 = it_tmp1 + 1;
        end
    end
    txt.XData = qa(:,1,i)+0.5;
    txt.YData = qa(:,2,i)+0.5;
%     for j = 1:Na
%         Gamma_Edge(j).XData = [qa(j,1,i) , qg(Aa(j),1,i)]; 
%         Gamma_Edge(j).YData = [qa(j,2,i) , qg(Aa(j),2,i)];
%     end
%     plot(time(i),my_dis(:,i),'k.','LineWidth',0.1);
    for j = 1:numel(d_L)
        circles(j).XData = qg(1,1,i) + d_L(j)*circle_x;
        circles(j).YData = qg(1,2,i) + d_L(j)*circle_y;       
    end  
    drawnow
%If you don't want to save output as GIF or EPS, you can comment the following commands       
    filename = strcat('obs_',num2str(i));
    print(filename,'-depsc2','-r600') %eps
    im = frame2im(getframe(f));
    [X,Map] = rgb2ind(im,256);
    imwrite(X , Map , 'obstacle.gif' , 'gif' , 'WriteMode','append' , 'Delay' , time_step);
end

%% Functions
%%% Equation 6
function y = sigmanorm(z , e) 
    y = (sqrt(1 + e*norm(z)^2)-1)/e;
end
%%% Equation 10
function y = a_ij(q , i , j , e , r , h) 
    y = rho(sigmanorm(q(j,:) - q(i,:) , e)/r , h);
end
%%% Equation 10
function y = b_ik(q , i , j , e , d , h) 
    y = rho(sigmanorm(q(j,:) - q(i,:) , e)/d , h);
end
%%% Equation 9 in "Flocking for Multi-Agent Dynamic Systems:
                %   Algorithms and Theory"
function y = sigmaepsilon(z , e) 
    y = z./sqrt(1 + e*(norm(z))^2);
end
%%% Equation 9
function y = n_ij(q , i , j , e) 
%     y = sigmaepsilon(q(j,:)-q(i,:) , e);
    y = (q(j,:)-q(i,:))/(sqrt(1+e*norm(q(j,:)-q(i,:))^2));
end
%%% Equation 9
function y = n_ri(q  , e) 
%     y = sigmaepsilon(q , e);
    y = (q /(sqrt(1+e*norm(q)^2)));
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
function y = phi_alpha(z , h , d , da , r)
    aa = 5;
    bb = 5;
    cc = abs(aa - bb)/sqrt(4*aa*bb);
    phi = @(z) 0.5*((aa+bb)*sigma_1(z/d+cc)+(aa-bb)); %%% Equation 7
    y = rho(z/r , h).*phi(z - da);
end
%%% Equation 7
function y = phi_gamma(z  ,d, d1 )
    aa = 5;
    bb = 5;
    cc = abs(aa - bb)/sqrt(4*aa*bb);
    phi = @(z) 0.5*((aa+bb)*sigma_1(z/d+cc)+(aa-bb)); %%% Equation 7
    y = 1*phi(z - d1);
end
%%% Equation 7
function y = phi_beta(z , h , db)
    y = rho(z/db , h).*(sigma_1(z - db)-1);
end
%%% Equation 8
function y = sigma_1(z)
    y = z/(1 + (z)^2)^0.5;
end
%%% Equation 4
function y = u(i , q , p , deltaq , deltap  , h  ...
    , d_L , d , da , da1_L , ra  ...
     , e  , c1a , c2a , c1b , c2b , c1g , c2g,yk, Rk, d_beta, q_hat, p_hat)
    y = [0,0];
    Na = size(q , 1);
    Nb = size(yk , 1);
    for j = 1:Na
        if j ~= i  
            y = y + c1a*phi_alpha(sigmanorm(q(j,:)-q(i,:),e) , h , d, da , ra ).*n_ij(q , i , j , e)...
                  + c2a*a_ij(q , i , j , e , ra , h)*(p(j,:)-p(i,:)); 
        end
    end
    for j = 1:Nb  
        y = y + c1b*phi_beta(sigmanorm(q_hat((Na*j-Na+i),:)-q(i,:),e) , h , d_beta).*n_ij([q ; q_hat] , i , (Na*j+i) , e)...
            + c2b*b_ik([q ; q_hat] , i , (Na*j+i) , e , d_beta , h)*(p_hat((Na*j-Na+i),:)-p(i,:));        
    end
    y = y - c1g*phi_gamma(sigmanorm(deltaq(i,:),e)  ,d_L , da1_L )... 
        .*n_ri(deltaq(i,:) , e) - c2g*deltap(i,:);
end
