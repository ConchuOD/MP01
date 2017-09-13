%% Preliminary Declarations & Calculations
% Resolution & Range of Theta0
Theta0 = [0:1:720]*(pi/180);
% Pre-allocations
Angles = zeros(length(Theta0), 10);
point0 = zeros(length(Theta0), 2);
pointa = zeros(length(Theta0), 2);
pointb = zeros(length(Theta0), 2);
pointc = zeros(length(Theta0), 2);
pointd = zeros(length(Theta0), 2);
pointe = zeros(length(Theta0), 2);
pointf = zeros(length(Theta0), 2);
pointg = ones(length(Theta0)+1, 2)*NaN;
M = struct('cdata', cell(1,length(Theta0)), 'colormap', cell(1,length(Theta0)));
% Defining lengths
a0 = 15;
ab = 50;
bd = 55.8;
de = 39.4;
eg = 65.7;
gf = 49;
ef = 36.7;
af = 61.9;
cf = 39.3;
cd = 40.1;
bc = 41.5;
% Initial Guesses
g1 = 150;
g2 = 210;
g3 = 80;
g4 = 180;
g5 = 270;
g6 = 180;
g7 = 290;
g8 = 250;
g9 = 260;
g10 = 300;
guesses = [g1 g2 g3 g4 g5 g6 g7 g8 g9 g10]*(pi/180);
% Counters
error_count = 0;
iteration_count = zeros(length(Theta0), 5);
% x and y coord of c
cx = -38;
cy = -7.8;
% Colours (normalised rgb)
barcolour = [222/256, 184/256, 135/256];
% Bar face matrix
facebar = [1 2 4 3];
% Setting up joint circles
ang = linspace(0,2*pi)';
rad = 1;
joint = [rad*cos(ang), rad*sin(ang)];
%% Calculate all angles,iteration_count(N,:),error_count
tic
for N = 1:length(Theta0)
    if N==1
    [Angles(N,:),error_count,iteration_count(N,:)] = Jansen(Theta0(N),guesses, error_count);
    else
    %
    [Angles(N,:),error_count,iteration_count(N,:)] = JansenV2(Theta0(N), Angles(N-1,:),error_count); 
    end
end
toc
%% From the angles, compute the positions of the points and draw each joint, which is a circle centred at the calculated point
tic
for N = 1:length(Theta0)
    % point0
    point0(N,:) = [0, 0];
    joint0(1,:,1) = joint(:,1);
    joint0(1,:,2) = joint(:,2); 
    % pointa
    pointa(N,:) = [a0*cos(Theta0(N)), a0*sin(Theta0(N))];
    jointa(1,:,1) = joint(:,1)+a0*cos(Theta0(N));
    jointa(1,:,2) = joint(:,2)+a0*sin(Theta0(N)); 
    % pointc
    pointc(N,:) = [cx, +cy];
    jointc(1,:,1) = joint(:,1)+cx;
    jointc(1,:,2) = joint(:,2)+cy;
    % pointb
    pointb(N,:) = [a0*cos(Theta0(N))+ab*cos(Angles(N,1)),a0*sin(Theta0(N))+ab*sin(Angles(N,1))];
    jointb(1,:,1) = joint(:,1)+a0*cos(Theta0(N))+ab*cos(Angles(N,1));
    jointb(1,:,2) = joint(:,2)+a0*sin(Theta0(N))+ab*sin(Angles(N,1));
    % pointd   
    pointd(N,:) = [cx+cd*cos(Angles(N,4)),cy+cd*sin(Angles(N,4))];
    jointd(1,:,1) = joint(:,1)+cx+cd*cos(Angles(N,4));
    jointd(1,:,2) = joint(:,2)+cy+cd*sin(Angles(N,4));
    % pointe
    pointe(N,:) = [cx+cd*cos(Angles(N,4))+de*cos(Angles(N,5)),+cy+cd*sin(Angles(N,4))+de*sin(Angles(N,5))];
    jointe(1,:,1) = joint(:,1)+cx+cd*cos(Angles(N,4))+de*cos(Angles(N,5));
    jointe(1,:,2) = joint(:,2)+cy+cd*sin(Angles(N,4))+de*sin(Angles(N,5));
    % pointf
    pointf(N,:) = [cx+cf*cos(Angles(N,7)),cy+cf*sin(Angles(N,7))];
    jointf(1,:,1) = joint(:,1)+cx+cf*cos(Angles(N,7));
    jointf(1,:,2) = joint(:,2)+cy+cf*sin(Angles(N,7));
    % pointg
    pointg(N,:) = [pointe(N,1)+eg*cos(Angles(N,10)),pointe(N,2)+eg*sin(Angles(N,10))];
    jointg(1,:,1) = joint(:,1)+pointe(N,1)+eg*cos(Angles(N,10));
    jointg(1,:,2) = joint(:,2)+pointe(N,2)+eg*sin(Angles(N,10));
%% Making bars
    % oa
    vertices0a = BarVertexGen(Theta0(N),point0(N,:),pointa(N,:)); 
    % ab
    verticesab = BarVertexGen(Angles(N,1),pointa(N,:),pointb(N,:));
    % bd
    verticesbd = BarVertexGen(Angles(N,2),pointd(N,:),pointb(N,:));
    % bc
    verticesbc = BarVertexGen(Angles(N,3),pointc(N,:),pointb(N,:));   
    % cd
    verticescd = BarVertexGen(Angles(N,4),pointc(N,:),pointd(N,:));
    % de
    verticesde = BarVertexGen(Angles(N,5),pointe(N,:),pointd(N,:));
    % ef
    verticesef = BarVertexGen(Angles(N,6),pointe(N,:),pointf(N,:));
    % cf
    verticescf = BarVertexGen(Angles(N,7),pointc(N,:),pointf(N,:));
    % fa
    verticesfa = BarVertexGen(Angles(N,8),pointa(N,:),pointf(N,:));
    % fg
    verticesfg = BarVertexGen(Angles(N,9),pointg(N,:),pointf(N,:));
    % eg
    verticeseg = BarVertexGen(Angles(N,10),pointe(N,:),pointg(N,:));       
%%  Visualisation  
    % Open a figure to plot the linkage
    jansenfig = figure('Position',[100 100 1600 900]);
    axis([-120 120 -100 100])
    % plot positions of foot
    patch('XData',pointg([1:N+1],1),'YData',pointg([1:N+1],2),'FaceColor','none','EdgeColor',[1 0 0]); hold on;
    % out of order on point 0 so crank is in background
    patch('Vertices',vertices0a,'Faces',facebar,'FaceColor',barcolour,'LineStyle','none'); hold on;
    patch(joint0(1,:,1),joint0(1,:,2),'k'); hold on;
    % patching bars
    patch('Vertices',verticesbd,'Faces',facebar,'FaceColor',barcolour,'LineStyle','none'); hold on;
    patch('Vertices',verticesab,'Faces',facebar,'FaceColor',barcolour,'LineStyle','none'); hold on;
    patch('Vertices',verticesbc,'Faces',facebar,'FaceColor',barcolour,'LineStyle','none'); hold on;
    patch('Vertices',verticescd,'Faces',facebar,'FaceColor',barcolour,'LineStyle','none'); hold on;
    patch('Vertices',verticesde,'Faces',facebar,'FaceColor',barcolour,'LineStyle','none'); hold on;
    patch('Vertices',verticesef,'Faces',facebar,'FaceColor',barcolour,'LineStyle','none'); hold on;
    patch('Vertices',verticescf,'Faces',facebar,'FaceColor',barcolour,'LineStyle','none'); hold on;
    patch('Vertices',verticesfg,'Faces',facebar,'FaceColor',barcolour,'LineStyle','none'); hold on;
    patch('Vertices',verticeseg,'Faces',facebar,'FaceColor',barcolour,'LineStyle','none'); hold on;
    patch('Vertices',verticesfa,'Faces',facebar,'FaceColor',barcolour,'LineStyle','none'); hold on;
    % patching joints
    patch(jointa(1,:,1),jointa(1,:,2),'k'); hold on;
    patch(jointb(1,:,1),jointb(1,:,2),'k'); hold on;
    patch(jointc(1,:,1),jointc(1,:,2),'k'); hold on;
    patch(jointd(1,:,1),jointd(1,:,2),'k'); hold on;
    patch(jointe(1,:,1),jointe(1,:,2),'k'); hold on;
    patch(jointf(1,:,1),jointf(1,:,2),'k'); hold on;
    patch(jointg(1,:,1),jointg(1,:,2),'r'); hold on;
    % Save frame to a movie
    axis off;
    M(N) = getframe(jansenfig);
    close(jansenfig)
end
toc
if(error_count ~= 0)
    warning('Newton Raphson did not converge %d times.', error_count)
end
y = [iteration_count([1:N],1) iteration_count([1:N],2) iteration_count([1:N],3) iteration_count([1:N],4) iteration_count([1:N],5)];
figure
bar(Theta0*(pi/180), y)
axis([ -0.5*(pi/180) 10.5*(pi/180) 0 7])
xlabel('Crank Angle (radians)')
ylabel('Number of iterations')
legend('Closure pair 1','Closure pair 2','Closure pair 4','Closure pair 3','Closure pair 5');

figure('Position',[100 100 1600 900]);
axis([-120 120 -100 100])
axis off
movie(M, 1)

% V = VideoWriter('Jansen_Animation_Final', 'MPEG-4');
% open(V)
% writeVideo(V,M)
% close(V)
V = VideoWriter('Jansen_Animation_14321611', 'Motion JPEG AVI');
open(V)
writeVideo(V,M)
close(V)
% V = VideoWriter('Jansen_Animation_14321611', 'Uncompressed AVI');
% open(V)
% writeVideo(V,M)
% close(V)
% clear all

