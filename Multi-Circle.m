clc
clear
close all
%% Initialize
% global my_phi
% simulation parameters
iterations = 2000;%200(sec)
time_step = 0.1;
tf = iterations*time_step; %final time
ts = time_step;
time = 0:ts:tf-0.1;
Tf=iterations;

k = 1.5;
d_L = 5;
e = 0.3;% Epsilon in Equation 5
h = 0.1;% h in Equation 10
%%% Parameters of Equation 4
c1a = 20;
c2a = 3;
c1g = 15;
c2g = 6;
%%%
N1 = 3; %Number of agent in first circle
N2 = 14; %Number of agent in second circle
N3 = 25; %Number of agent in third circle
N4 = 36; %Number of agent in fourth circle
d = (sin(pi/(N1))*d_L*2);%Equation (14)  Desired Distance Between Agents
d_L2= (sin(pi/(N1))*d_L*2)/(2*sin(pi/N2));% Desired Radius of Second Circle
d_L3= (sin(pi/(N1))*d_L*2)/(2*sin(pi/N3));% Desired Radius of Third Circle
d_L4= (sin(pi/(N1))*d_L*2)/(2*sin(pi/N4));% Desired Radius of Fourth Circle
dL = [d_L;d_L2;d_L3;d_L4];


r = d*k
ra   = sigmanorm(r,e);
da   = sigmanorm(d,e);


kk = 1;
r1_L= d_L+d_L/kk;%Attraction Range (Cut-Off) of First Circle 
r2_L= d_L2+d_L/kk;
r3_L= d_L3+d_L/kk;
r4_L= d_L4+d_L/kk;
rL = [r1_L;r2_L;r3_L;r4_L];

ra1_L = sigmanorm(r1_L,e);
ra2_L = sigmanorm(r2_L,e);
ra3_L = sigmanorm(r3_L,e);
ra4_L = sigmanorm(r4_L,e);

da1_L = sigmanorm(d_L,e);
da2_L = sigmanorm(d_L2,e);
da3_L = sigmanorm(d_L3,e);
da4_L = sigmanorm(d_L4,e);


% if (d_L2-r1_L)<r||(d_L3-r2_L)<r||(d_L4-r3_L)<r
%     disp(['Error: Interaction range between agents is greater than one of distance between two circles'])
%     return
% end
epsilin1 = 0.1;
d_L2-d_L-epsilin1
d_L3-d_L2-epsilin1
d_L4-d_L3-epsilin1
if (d_L2-d_L-epsilin1)<r||(d_L3-d_L2-epsilin1)<r||(d_L4-d_L3-epsilin1)<r
    disp(['Error: Interaction range between agents is greater than one of distance between two circles'])
    return
end
% define gamma agents

Ng = 1;

qg = 0*ones(Ng , 2 , iterations);
pg = zeros(Ng , 2 , iterations);

pg(1,:,:) = 2*zeros(1 , 2 , iterations);
start_L = 10;
pg(1,:,start_L:end) = 2*ones(1 , 2 , iterations-start_L+1);
% define alpha agents

Na = N1+N2+N3+N4;
qa = zeros(Na , 2 , iterations);
pa = zeros(Na , 2 , iterations);

%%
%Initial Position of Agents
% Generate unique coordinates between the circles

generatedPoints = zeros(N1, 2);
count = 0;
EP = 0.0;
minDistance = 4;
while count < N1
%     r = randi([0, d_L]);  % Random integer radius
    r = (r1_L - 0) * rand() + 0;  % Random non-integer radius
    theta = rand() * 2 * pi;  % Random angle between 0 and 2*pi
    x_temp = r * cos(theta);  %
    y_temp = r * sin(theta);  
    
    % Check if the point is inside the outer circle but outside the inner circle

    if r > 0 && r < (r1_L + EP)
        % Check if the point is unique
        if ~any(ismember(generatedPoints, [x_temp, y_temp], 'rows'))&& ...
            all(sqrt((qa(:,1,1) - x_temp).^2 + (qa(:,2,1) - y_temp).^2) >= 0.1)
            count = count + 1;          
            qa(count,1,1) = x_temp;
            qa(count,2,1) = y_temp;
            generatedPoints(count, :) = [x_temp, y_temp];
        end
    end
end
while count < N1+N2
%     r = randi([d_L, d_L2]);  % Random integer radius
    r = (r2_L - r1_L) * rand() + r1_L;  % Random non-integer radius
    theta = rand() * 2 * pi;  % Random angle between 0 and 2*pi
%     x_temp = round(r * cos(theta));  % Integer x-coordinate
%     y_temp = round(r * sin(theta));  % Integer y-coordinate
    x_temp = r * cos(theta);  
    y_temp = r * sin(theta);     
    % Check if the point is inside the outer circle but outside the inner circle
    if r > (r1_L + EP) && r < (r2_L + EP)
        % Check if the point is unique
        if ~any(ismember(generatedPoints, [x_temp, y_temp], 'rows'))&& ...
            all(sqrt((qa(:,1,1) - x_temp).^2 + (qa(:,2,1) - y_temp).^2) >= minDistance)
            count = count + 1;
            qa(count,1,1) = x_temp;
            qa(count,2,1) = y_temp;
            generatedPoints(count, :) = [x_temp, y_temp];
        end
    end
end
while count < N1+N2+N3
%     r = randi([d_L2, d_L3]);  % Random integer radius
    r = (r3_L - r2_L) * rand() + r2_L;  % Random non-integer radius
    theta = rand() * 2 * pi;  % Random angle between 0 and 2*pi
%     x_temp = round(r * cos(theta));  % Integer x-coordinate
%     y_temp = round(r * sin(theta));  % Integer y-coordinate
    x_temp = r * cos(theta); 
    y_temp = r * sin(theta);
    % Check if the point is inside the outer circle but outside the inner circle
    if r > (r2_L + EP) && r < (r3_L + EP)
        % Check if the point is unique
        if ~any(ismember(generatedPoints, [x_temp, y_temp], 'rows'))&& ...
            all(sqrt((qa(:,1,1) - x_temp).^2 + (qa(:,2,1) - y_temp).^2) >= minDistance)
            count = count + 1;
            qa(count,1,1) = x_temp;
            qa(count,2,1) = y_temp;
            generatedPoints(count, :) = [x_temp, y_temp];
        end
    end
end
while count < N1+N2+N3+N4
%     r = randi([d_L3, d_L4]);  % Random integer radius
    r = (r4_L - r3_L) * rand() + r3_L;  % Random non-integer radius
    theta = rand() * 2 * pi;  % Random angle between 0 and 2*pi
%     x_temp = round(r * cos(theta));  % Integer x-coordinate
%     y_temp = round(r * sin(theta));  % Integer y-coordinate
    x_temp = r * cos(theta);
    y_temp = r * sin(theta);  
    % Check if the point is inside the outer circle but outside the inner circle
    if r > (r3_L + EP) %&& r < r4_L
        % Check if the point is unique
        if ~any(ismember(generatedPoints, [x_temp, y_temp], 'rows'))&& ...
            all(sqrt((qa(:,1,1) - x_temp).^2 + (qa(:,2,1) - y_temp).^2) >= minDistance)
            count = count + 1;
            qa(count,1,1) = x_temp;
            qa(count,2,1) = y_temp;
            generatedPoints(count, :) = [x_temp, y_temp];
        end
    end
end

%%
N_link = 0.5*Na*(Na-1);
distance_qa = zeros(N_link , iterations);
Dis_Q = [];
distance_AG = zeros(Na,iterations); % Leader & Agents
for i = 1:Na
distance_AG(i,1) = i;
end
ErrorAG = zeros(Na,iterations);
uc = zeros(Na , 2 , iterations);
%% simulation

deltaq = zeros(Na , 2);
deltap = zeros(Na , 2);

it = 1;
while it <= iterations
 
    % calculate qdelta and pdelta
    QG = qg(:,:,it);
    PG = pg(:,:,it);
    QA = qa(:,:,it);
    PA = pa(:,:,it);
    for i = 1:Na
        deltaq(i,:) = QA(i,:) - QG;
        deltap(i,:) = PA(i,:) - PG;
    end
    for i = 1:Na        
            uc(i,:,it) = u(i , QA , PA ,...
                deltaq , deltap ,...
                h , d_L ,  da , da1_L, da2_L ,da3_L ,da4_L ...
                , ra , ra1_L , ra2_L ,ra3_L , e ,...
                c1a , c2a , c1g , c2g  , d );      
    end
    pa(:,:,it+1) = pa(:,:,it) + time_step*(uc(:,:,it));
    qg(:,:,it+1) = qg(:,:,it) + time_step*pg(:,:,it);
%     pg(:,:,it+1) = pg(:,:,it);
    qa(:,:,it+1) = qa(:,:,it) + time_step*pa(:,:,it);
    %%% Total Distance Calculate
%     number_link = 1;
% %     Dis_Q = distance_qa(:,it);
%     it1 = it+1;
%     Dis_Q = [];
%     QGJ = qg(1,:,it1);
%     for i=1:Na
%         QAi = QA(i,:);    
%         DIS = zeros(Na-i,1);
%         for j=(i+1):Na
%                 QAJ = qa(j,:,it);
%                 DIS(j-i) = norm(QAi-QAJ);
%         end
% %         distance_qa(Na*(i-1)-(0.5*(i-2)*(i+1)):i*Na-0.5*i(i+1),it) = DIS;
%         Dis_Q = [Dis_Q;DIS];
%         
%         distance_AG(i,it1) = norm(qa(i,:,it1)-QGJ);
%     end
%     distance_qa(:,it) = Dis_Q;

    %%% Total Distance Calculate Method2
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
%     Min_dis(:,it) = sortrows(distance_qa(:,it));


    %%%
    it = it + 1;
end
%%% Distance between Leader and Followers
distance_AG = sortrows(distance_AG,iterations);
%%% Distance Between Agents 
my_dis = sortrows(distance_qa,iterations); % Distance Between Agents in First Circle

my_dis = my_dis(1:Na,:);


%% Final Plot
figure;
% subplot(3,1,1)
circle_x = cosd(1:360);
circle_y = sind(1:360);

Agents = scatter(qa(:,1,Tf) , qa(:,2,Tf) , 'bo');
hold on
Leaders = scatter(qg(:,1,Tf) , qg(:,2,Tf) , 'r*');


Alpha_Edge = [];
for i = 1:Na
    for j = i+1:Na
        Alpha_Edge =  plot([qa(i,1,Tf) , qa(j,1,Tf)] , [qa(i,2,Tf) , qa(j,2,Tf)]); 
        if a(qa(:,:,Tf) , i , j , e , ra , h) == 0
            Alpha_Edge(end).Color = 'none';
        else
            Alpha_Edge(end).Color = 'm';
        end
    end
end
circles2 = [];
for i = 1:numel(rL)-1
    circles2 = [circles2 plot(qg(1,1,Tf) + rL(i)*circle_x , qg(1,2,Tf) + rL(i)*circle_y , 'g-')];
end
% circles1 = [];
% for i = 1:numel(dL)
%     circles1 = [circles1 plot(qg(1,1,Tf) + dL(i)*circle_x , qg(1,2,Tf) + dL(i)*circle_y , 'k-')];
% end
% str_a = string(1:Na);
% txt = textscatter(2+qa(:,1,Tf),2+qa(:,2,Tf),str_a);


T = title('$t = 0.00 $(sec) ' ,'interpreter','latex','FontSize',14);
T.String = sprintf("t = %2.2f (sec) " , time_step*(Tf));
xlabel('$X$' ,'interpreter','latex','FontSize',14)
ylabel('$Y$' ,'interpreter','latex','FontSize',14)
axis equal

figure;
plot(time,distance_AG);
title('All Distances Between Leader \& Agents','interpreter','latex','FontSize',14)
xlabel('Time $(s)$','interpreter','latex','FontSize',14)
ylabel('Distance $(m)$','interpreter','latex','FontSize',14)
axis ([0 tf ,0 60])
print('AL_M','-depsc2','-r600')
figure;
plot(time,my_dis);
title('Shortest Distances Between Agents','interpreter','latex','FontSize',14)
xlabel('Time $(s)$','interpreter','latex','FontSize',14)
ylabel('Distance $(m)$','interpreter','latex','FontSize',14)
axis ([0 tf ,0 16])
print('AA_M','-depsc2','-r600')

ErrorAG(1:N1,:) = distance_AG(1:N1,:)-d_L;
ErrorAG(N1+1:N2+N1,:) = distance_AG(N1+1:N2+N1,:)-d_L2;
ErrorAG(N2+N1+1:N2+N1+N3,:) = distance_AG(N2+N1+1:N2+N1+N3,:)-d_L3;
ErrorAG(N2+N1+N3+1:N2+N1+N3+N4,:) = distance_AG(N2+N1+N3+1:N2+N1+N3+N4,:)-d_L4;
% 
figure;
Leader_Error = rms(ErrorAG);
P1 = plot(time,Leader_Error);
title('RMSE of the Distances Between Leader \& Agents','interpreter','latex','FontSize',14);
xlabel('Time $(s)$','interpreter','latex','FontSize',14)
ylabel('Value $(m)$','interpreter','latex','FontSize',14)
saveas(P1, 'AL_Error_M.png')
print('AL_Error_M','-depsc2','-r600')
% % axis ([0 tf ,0 20])
% 
Adj_Error = rms(my_dis-d);
figure;
P2 = plot(time,Adj_Error);
title('RMSE of the Shortest Distances Between Agents','interpreter','latex','FontSize',14);
xlabel('Time $(s)$','interpreter','latex','FontSize',14)
ylabel('Value $(m)$','interpreter','latex','FontSize',14)
saveas(P2, 'AA_Error_M.png')
print('AA_Error_M','-depsc2','-r600')
% 
steady_error = iterations - 799;
figure;
P3 = plot(time(:,steady_error:end),Adj_Error(:,steady_error:end));
saveas(P3, 'AA_Zoom_M.png')
print('AA_Zoom_M','-depsc2','-r600')
figure;
P4 = plot(time(:,steady_error:end),Leader_Error(:,steady_error:end));
print('AL_Zoom_M','-depsc2','-r600')
saveas(P4, 'AL_Zoom_M.png')

% axis ([0 tf ,0 20])

% % 
% % tss = 700;
% % nf = iterations - tss+1;
% % my_dis_f2 = my_dis(1:Na,tss:end);
% % D_L1 = distance_AG(1:N1,tss:end);
% % D_L2 = distance_AG(N1+1:N2+N1,tss:end);
% % D_L3 = distance_AG(N2+N1+1:N2+N1+N3,tss:end);
% % D_L4 = distance_AG(N2+N1+N3+1:N2+N1+N3+N4,tss:end);
% % Actual_DL = [D_L1;D_L2;D_L3;D_L4];
% % Desired_DL = [d_L.*ones(N1,nf);d_L2.*ones(N2,nf);d_L3.*ones(N3,nf);d_L4.*ones(N4,nf)];
% % meanDL = d_L+d_L2+d_L3+d_L4/numel(dL);
% % X = [(meanDL/d)*my_dis_f2 Actual_DL];
% % XD = [(meanDL/d)*d.*ones(Na,nf) Desired_DL];
% % F = norm(XD - X,'fro')
% % end
f = figure;
% subplot(2,1,1)

circle_x = cosd(1:360);
circle_y = sind(1:360);

Agents = scatter(qa(:,1,1) , qa(:,2,1) , 'bo');
hold on
Leaders = scatter(qg(:,1,1) , qg(:,2,1) , 'r*');

Alpha_Edge = [];
for i = 1:Na
    for j = i+1:Na
        Alpha_Edge = [Alpha_Edge plot([qa(i,1,1) , qa(j,1,1)] , [qa(i,2,1) , qa(j,2,1)])]; 
        if a(qa(:,:,1) , i , j , e , ra , h) == 0
            Alpha_Edge(end).Color = 'none';
        else
            Alpha_Edge(end).Color = 'm';
        end
    end
end

% Gamma_Edge = [];
% for i = 1:Na
%     Gamma_Edge = [Gamma_Edge plot([qa(i,1,1) , qg(Aa(i),1,1)] , [qa(i,2,1) , qg(Aa(i),2,1)] , 'g-')]; 
% end

%%%
circles = [];
for i = 1:numel(dL)
    circles = [circles plot(qg(1,1,1) + dL(i)*circle_x , qg(1,2,1) + dL(i)*circle_y , 'k--','LineWidth',0.25)];
end
circles2 = [];
for i = 1:numel(rL)
    circles2 = [circles2 plot(qg(1,1,1) + rL(i)*circle_x , qg(1,2,1) + rL(i)*circle_y , 'g-')];
end
%%%

T = title('t = 0.00 (sec)','interpreter','latex','FontSize',14);
xlabel('$x$','interpreter','latex','FontSize',14)
ylabel('$y$','interpreter','latex','FontSize',14)
axis equal

im = frame2im(getframe(f));
[X,Map] = rgb2ind(im,256);
imwrite(X , Map , 'muli-circle-formation.gif' , 'gif' , 'Loopcount',inf , 'Delay' , time_step);
%%
% zm = d_L*4;
for i = 1:iterations
    if i ==1
%     circles = [];
%     circles2= [];
        print('multi_3_init','-depsc2','-r600')
    end  
    T.String = sprintf("t = %2.2f (sec)" , time_step*(i));
    Agents.XData = qa(:,1,i+1);
    Agents.YData = qa(:,2,i+1);
    
    Leaders.XData = qg(:,1,i+1);
    Leaders.YData = qg(:,2,i+1);
    
    for j = 1:numel(rL)
        circles(j).XData = qg(1,1,i+1) + dL(j)*circle_x;
        circles(j).YData = qg(1,2,i+1) + dL(j)*circle_y;
        circles2(j).XData = qg(1,1,i+1) + rL(j)*circle_x;
        circles2(j).YData = qg(1,2,i+1) + rL(j)*circle_y;        
    end    
    it_tmp = 1;
    for j = 1:Na
        for k = j+1:Na
            Alpha_Edge(it_tmp).XData = [qa(j,1,i+1) , qa(k,1,i+1)];
            Alpha_Edge(it_tmp).YData = [qa(j,2,i+1) , qa(k,2,i+1)];
            if a(qa(:,:,i) , j , k , e , ra , h) == 0
                Alpha_Edge(it_tmp).Color = 'none';
            else
                Alpha_Edge(it_tmp).Color = 'm';
            end
            
            it_tmp = it_tmp + 1;
        end
    end
    
%     for j = 1:Na
%         Gamma_Edge(j).XData = [qa(j,1,i) , qg(Aa(j),1,i)];
%         Gamma_Edge(j).YData = [qa(j,2,i) , qg(Aa(j),2,i)];
%     end
if i ==iterations
%     circles = [];
%     circles2= [];
    print('multi_3_final','-depsc2','-r600')
end

    drawnow
%If you don't want to save output as GIF or EPS, you can comment the following commands       

    im = frame2im(getframe(f));
    [X,Map] = rgb2ind(im,256);
    imwrite(X , Map , 'muli-circle-formation.gif' , 'gif' , 'WriteMode','append' , 'Delay' , time_step);
end
%% Functions
%%% Equation 6
function y = sigmanorm(z , e) 
    y = (sqrt(1 + e*norm(z)^2)-1)/e;
end
%%% Equation 10
function y = a(q , i , j , e , r , h) 
    y = rho(sigmanorm(q(j,:) - q(i,:) , e)/r , h);
end
function y = sigmaepsilon(z , e) 
    y = z./sqrt(1 + e*(norm(z))^2);
end
%%% Equation 9
function y = n_ij(q , i , j , e) 
    y = sigmaepsilon(q(j,:)-q(i,:) , e);
end
%%% Equation 9
function y = n_ri(q , e) 
    y = sigmaepsilon(q , e);
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
function y = phi_alpha(z , h , d , r, D)
    aa = 5;
    bb = 5;
    cc = abs(aa - bb)/sqrt(4*aa*bb);
    phi = @(z) 0.5*((aa+bb)*sigma_11(z/D+cc)+(aa-bb));
    y = rho(z/r , h).*phi(z - d);
end
%%% Equation 15
function y = phi_alpha_gamma(z , d , d1 , d2 , d3 , d4 , r1 , r2 , r3)
    aa = 5;
    bb = 5;
    cc = abs(aa - bb)/sqrt(4*aa*bb);
    phi1 = @(z) 0.5*((aa+bb)*sigma_11(z/d+cc)+(aa-bb));
        if     z<=r1
            y = 1*phi1(z - d1);
        elseif z<=r2
            y = 1*phi1(z - d2);
        elseif z<=r3
            y = 1*phi1(z - d3);
        else
            y = 1*phi1(z - d4);
        end
end
%%% Equation 8
function y = sigma_11(z)
    y = z/(1 + (z)^2)^0.5;
end
%%% Equation 4
function y = u(i , q , p , deltaq , deltap , h  , d_L , da , da1_L, da2_L , da3_L , da4_L , ra , ra1_L , ra2_L , ra3_L , e  , c1a , c2a , c1g , c2g, d)
    y = [0,0];
    Na = size(q , 1);

    for j = 1:Na
        if j ~= i
            y = y + c1a*phi_alpha(sigmanorm(q(j,:)-q(i,:),e) , h , da , ra , d).*n_ij(q , i , j , e)...
                  + c2a*a(q , i , j , e , ra , h)*(p(j,:)-p(i,:)); 
        end
    end

    y = y - c1g*phi_alpha_gamma(sigmanorm(deltaq(i,:),e) , d_L, da1_L , da2_L , da3_L , da4_L, ra1_L , ra2_L, ra3_L)... 
        .*n_ri(deltaq(i,:) , e) - c2g*deltap(i,:);
end