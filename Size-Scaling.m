close all
clear
clc
%% Initialize
% global my_phi
% simulation parameters

iterations = 1200;
time_step = 0.1;
tf = iterations*time_step; %final time
ts = time_step;
time = 0:ts:tf-0.1;
Tf=iterations;
%%%
Na = 12;
k = 1.3; 
e = 0.3; % Epsilon in Equation 6
h = 0.1; % h in Equation 11
d_L = 5; %Radius of circle formation
d_gamma = sigmanorm(d_L,e);
d = sin(pi/Na)*d_L*2;% Equation (14),  Desired Distance Between Agents
d_alpha=sigmanorm(d,e);
r = k*d;
r_alpha= sigmanorm(r,e);

%%% Parameters of Equation 4
c1a = 10;
c2a = 4;
c1b = 40;
c2b = 10;
c1g = 7;
c2g = 5;
% define gamma agents

Ng = 1;% Number of leader

qg = 0*ones(Ng , 2 , iterations);
pg = zeros(Ng , 2 , iterations);
start_L = 400;
pg(1,:,:) = 2*zeros(1 , 2 , iterations);
pg(1,:,start_L:end) = 2*ones(1 , 2 , iterations-start_L+1);

% define alpha agents

qa = zeros(Na , 2 , iterations);
% qa(:,:,1) = [randperm(10,Na) ;randperm(10,Na)]';
qa(:,:,1) = 10*randn(Na , 2 );



pa = zeros(Na , 2 , iterations);
Aa = randi(Ng , Na , 1);

N_link = 0.5*Na*(Na-1);
distance_qa = zeros(N_link , iterations);
distance_AG = zeros(Na,iterations); % Leader & Agents
for i = 1:Na
distance_AG(i,1) = i;
end
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

uc = zeros(Na , 2 , iterations);

test_obs = zeros(Na^2 , 2 , iterations);
% test1_obs = zeros(Na^2-Na , 2 , iterations+1);
test2_obs= zeros(Na*Nb , 2 , iterations);

obs1 = zeros(Na , 2 , iterations);
obs2 = zeros(2*Na , 2 , iterations);
%% simulation

deltaq = zeros(Na , 2);% qa - qg
deltap = zeros(Na , 2);
it = 1;
num_flt = Na;
time_flt = 300;
    %%%
    
%%%
adj2 =[];
wait = 0;
wait1= 0;
wait2 = 0;
my_dis1 = zeros(N_link,iterations);
dis_AG = zeros(Na,iterations);
DA = zeros(1,iterations);
DL = zeros(1,iterations);
while it <= iterations

    for i = 1:Na
        deltaq(i,:) = qa(i,:,it) - qg(Aa(i),:,it);
        deltap(i,:) = pa(i,:,it) - pg(Aa(i),:,it);
    end
    % update
    qg(:,:,it+1) = qg(:,:,it) + time_step*pg(:,:,it);
    qa(:,:,it+1) = qa(:,:,it) + time_step*pa(:,:,it);


    number_link = 1;
    number_link1 = 1;
    number_link2 = 1;
    for i = 1:Na
        for j = 1:Nb
            %%% Equation 5
            mu = Rk(j)/(norm(qa(i,:,it)-yk(j,:,it)));
            ak = (qa(i,:,it)-yk(j,:,it))/(norm(qa(i,:,it)-yk(j,:,it)));
            P = eye - ak*transpose(ak);
            q_hat((Na*j-Na+i),:,it) = mu*qa(i,:,it)+(1-mu)*yk(j,:,it);
            p_hat((Na*j-Na+i),:,it) = mu*P*pa(i,:,it);
            %%%
        end
        uc(i,:,it) = u(i , qa(:,:,it) , pa(:,:,it) ,...
                                       deltaq , deltap ,...
                                       h  , d_L , d , d_alpha , d_gamma ...
                                       , r_alpha   , e  , c1a , c2a , c1b , c2b , c1g , ...
                                       c2g, yk(:,:,it), Rk, d_beta, ...
                                       q_hat(:,:,it),p_hat(:,:,it));     
    end
    pa(:,:,it+1) = pa(:,:,it) + time_step*(uc(:,:,it));% + pg(Aa(:),:,it+1) - pg(Aa(:),:,it));
    %%% Total Distance Calculate
    for i = 1:Na
        for j = 1:Na
            if j ~= i
                test_obs(number_link1,:,it)=a_ij(qa(:,:,it) , i , j , e , r_alpha , h);               
                number_link1 = number_link1 +1;
            else
                test_obs(number_link1,:,it)=0;
                number_link1 = number_link1 +1;
            end
        end
        for j = 1:Nb
            test2_obs(number_link2,:,it)=b_ik([qa(:,:,it) ; q_hat(:,:,it)] , i , (Na*j+i) , e , d_beta , h);
            number_link2 = number_link2+1;
        end
        for j=i+1:Na
            distance_qa(number_link,it) = norm(qa(i,:,it)-qa(j,:,it));
            number_link = number_link+1;
        end
        distance_AG(i,it) = norm(qa(i,:,it)-qg(1,:,it));
    end
    my_dis1(:,it) = sortrows(distance_qa(:,it));
    dis_AG(:,it) = sortrows(distance_AG(:,it));
    %%%
    %%%
    if it>start_L
    for i=1:Na
        obs1(i,:,it) = test2_obs(3*i-2,:,it);
        obs2(2*i-1:2*i,:,it) = test2_obs(3*i-1:3*i,:,it);
    end
    end
    adj2 =[];
    obs  = obs1(:,1,it);
    obs_ = obs2(:,1,it);

    if (numel(obs(obs>0))>=(Na/4)) %&& (wait ==0 )
        wait = wait+1;
%         disp('ok')
        indices = find(obs);
        for i=1:numel(indices)
        adj = find(test_obs(Na*indices(i)-(Na-1):Na*indices(i),1,it));
        adj2 = [adj2;sum(ismember(indices,adj))];
        end
    end    

    if all(adj2 ~=0) && (wait >= 3 )
%             disp('ok')
        d_L1 = 10;
        d_L  = d_L1;
        d_gamma = sigmanorm(d_L,e);
        time_r = it;
        wait = 0;
        wait1= 100;
        d = sin(pi/Na)*d_L*2;
        d_alpha=sigmanorm(d,e);
        r = k*d;
        r_alpha= sigmanorm(r,e);
    end

    if (numel(obs_(obs_>0))>=(Na/4)) %&& (wait2 ==0 )
%         disp('ok')
        indices1 = find(obs_(1:Na));
        indices2 = find(obs_(Na+1:end));
        indices = [indices1;indices2];
        for i=1:numel(indices)
        adj = find(test_obs(Na*indices(i)-(Na-1):Na*indices(i),1,it));
        adj2 = [adj2;sum(ismember(indices,adj))];
        end
        wait2 = wait2 + 1;
    end  

    if ((wait2>=3)  && (numel(adj2) ~= nnz(adj2)))
        d_L2= 2.5;
        d_L = d_L2;
        d_gamma = sigmanorm(d_L,e);
        wait2 = 0;
        time_r = [time_r;it];
        wait1 = 50;
        d = sin(pi/Na)*d_L*2;
        d_alpha=sigmanorm(d,e);
        r = k*d;
        r_alpha= sigmanorm(r,e);
        
    end

%     if (numel(obs(obs>0)))== 0
%                     d_L = 5;
%             da1_L = sigmanorm(d_L,e_L);
%             d = sin(pi/Na)*d_L*2;
%             da=sigmanorm(d,e);
%             r = k*d;
%             ra= sigmanorm(r,e);
%             time_r = [time_r;it];
%     end
     if wait1 ~=0      
        wait1=wait1-1;
        wait = 0;
        wait2 = 0;
        if wait1 == 0
            d_L = 5;
            d_gamma = sigmanorm(d_L,e);
            d = sin(pi/Na)*d_L*2;
            d_alpha=sigmanorm(d,e);
            r = k*d;
            r_alpha= sigmanorm(r,e);
            time_r = [time_r;it];
            
        end
    end
    %%%
    DA(1,it) = d;
    DL(1,it) = d_L;
    it = it + 1;    
end

dL = [d_L1 d_L d_L2 d_L];
%% plot Distance
figure;
plot(time,distance_qa);
title('All Distance Between Each Agent')
xlabel('Time')
ylabel('Distance')
axis ([0 tf ,0 50])

figure;
plot(time,dis_AG);
title('All Distances Between Leader \& Agents','interpreter','latex','FontSize',14)
xlabel('Time $(s)$','interpreter','latex','FontSize',14)
ylabel('Distance  $(m)$','interpreter','latex','FontSize',14)
axis ([0 tf ,0 20])
print('AL_Scl','-depsc2','-r600')

figure;
plot(time,my_dis1(1:Na,:));
title('Shortest Distances Between Agents','interpreter','latex','FontSize',14)
xlabel('Time $(s)$','interpreter','latex','FontSize',14)
ylabel('Distance $(m)$','interpreter','latex','FontSize',14)
axis ([0 tf ,0 10])
print('AA_Scl','-depsc2','-r600')
%%% Error Calculation
AL_Error = rms(DL - dis_AG);
AA_Error = rms(DA - my_dis1(1:Na,:));
starter_error = iterations - 199;
figure;
P1 = plot(time,AL_Error);
title('RMSE of the Distances Between Leader \& Agents','interpreter','latex','FontSize',14)
xlabel('Time $(s)$','interpreter','latex','FontSize',14)
ylabel('Value $(m)$','interpreter','latex','FontSize',14)
axis ([0 tf ,0 15])
print('AL_Error_Scl','-depsc2','-r600')
saveas(P1, 'AL_Error_Scl.png')
figure;
P2 = plot(time(1,starter_error:end),AL_Error(:,starter_error:end));
print('AL_Zoom_Scl','-depsc2','-r600')
saveas(P2, 'AL_Zoom_Scl.png')
figure;
P3 = plot(time,AA_Error);
title('RMSE of the Shortest Distances Between Agents','interpreter','latex','FontSize',14)
xlabel('Time $(s)$','interpreter','latex','FontSize',14)
ylabel('Value $(m)$','interpreter','latex','FontSize',14)
axis ([0 tf ,0 10])
print('AA_Error_Scl','-depsc2','-r600')
saveas(P3, 'AA_Error_Scl.png')
figure;
P4 = plot(time(1,starter_error:end),AA_Error(:,starter_error:end));
print('AA_Zoom_Scl','-depsc2','-r600')
saveas(P4, 'AA_Zoom_Scl.png')
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

% Gamma_Edge = [];
% for i = 1:Na
%     Gamma_Edge = [Gamma_Edge plot([qa(i,1,Tf) , qg(Aa(i),1,Tf)] , [qa(i,2,Tf) , qg(Aa(i),2,Tf)] , 'g-')]; 
% end

T = title('t = 0.00 (sec)');
T.String = sprintf("t = %2.2f (sec) " , time_step*(Tf) );
xlabel('X')
ylabel('Y')
axis equal

% subplot(3,1,2)
% plot(time,distance_AG);
% title('All Distance Between Leader & Agents')
% xlabel('Time')
% ylabel('Distance')
% axis ([0 tf ,0 8])
%     grid on 
%     grid minor
% subplot(3,1,3)
% if  my_dis>=0
%     plot(time,my_dis);
%     title('My Desired Distance')
%     xlabel('Time')
%     ylabel('Distance')
%     grid on 
%     grid minor
%     axis ([0 tf ,0 3])
% end

%% Animation Plot

f = figure;
% subplot(2,1,1)
        
        d = sin(pi/12)*5*2;
        d_alpha=sigmanorm(d,e);
        r = 1.2*d;
        r_alpha= sigmanorm(r,e);
        
circle_x = cosd(1:360);
circle_y = sind(1:360);

Agents = scatter(qa(:,1,1) , qa(:,2,1) , 'bo');
hold on
Leaders = scatter(qg(:,1,1) , qg(:,2,1) , 'r*');
Obstacles = [];
for i = 1:Nb
    Obstacles = [Obstacles plot(yk(i,1,1) + Rk(i)*circle_x , yk(i,2,1) + Rk(i)*circle_y , 'k-')];
end

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
% for i = 1:Na
%     Gamma_Edge = [Gamma_Edge plot([qa(i,1,1) , qg(Aa(i),1,1)] , [qa(i,2,1) , qg(Aa(i),2,1)] , 'g-')]; 
% end
circles = [];

circles = [circles plot(qg(1,1,1) + dL(1)*circle_x , qg(1,2,1) + dL(1)*circle_y , 'k--','LineWidth',0.25)];

T = title('$t = 0.00$ (sec)','interpreter','latex','FontSize',14);
xlabel('$x$','interpreter','latex','FontSize',14)
ylabel('$y$','interpreter','latex','FontSize',14)
axis equal


% subplot(2,1,2)
% hold on
% title('Distance between Agents')
% xlabel('Time')
% ylabel('Distance Vakue')
% dis = [];
% for j = 1:size_r
%     dis = [dis plot(time(1),my_dis(j,1),'k.','LineWidth',0.1)];
% end
% grid on 
% grid minor
% axis([0,30,0,60]);
% axis equal


im = frame2im(getframe(f));
[X,Map] = rgb2ind(im,256);
imwrite(X , Map , 'Scaling.gif' , 'gif' , 'Loopcount',inf , 'Delay' , time_step);
%%
tester1=1;

for i = 1:iterations
    if i ==1
%     circles = [];
%     circles2= [];
        print('scaling_0','-depsc2','-r600')
    end
    if i>(start_L+5)
        axis([qg(:,1,i)-(d_L+8) qg(:,1,i)+(d_L+8),qg(:,2,i)-(d_L+8) qg(:,2,i)+(d_L+8)])
    end
%     if any(ismember(time_r,i-2))
%         axis([qg(:,1,i)-(d_L+5) qg(:,1,i)+(d_L+5),qg(:,2,i)-(d_L+5) qg(:,2,i)+(d_L+5)])
%     end
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
    for j=1:numel(time_r)
        if i == time_r(j) + 2
            d = sin(pi/Na)*dL(j)*2;
            d_L = dL(j);
            d_alpha=sigmanorm(d,e);
            r = 1.2*d;
            r_alpha= sigmanorm(r,e);
        end
    end
    it_tmp1 = 1;
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
%     for j = 1:Na
%         Gamma_Edge(j).XData = [qa(j,1,i) , qg(Aa(j),1,i)]; 
%         Gamma_Edge(j).YData = [qa(j,2,i) , qg(Aa(j),2,i)];
%     end
%     plot(time(i),my_dis(:,i),'k.','LineWidth',0.1);

    circles.XData = qg(1,1,i) + d_L*circle_x;
    circles.YData = qg(1,2,i) + d_L*circle_y;       
 
    drawnow

    filename = strcat('scaling_',num2str(i));
    print(filename,'-depsc2','-r600')
    im = frame2im(getframe(f));
    [X,Map] = rgb2ind(im,256);
    imwrite(X , Map , 'Scaling.gif' , 'gif' , 'WriteMode','append' , 'Delay' , time_step);
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
