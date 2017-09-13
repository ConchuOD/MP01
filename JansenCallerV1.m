% Resolution
Theta0 = [0:10:720]*(pi/180);
% Pre-allocation of matrices
Angles = zeros(length(Theta0), 10);
point0 = zeros(length(Theta0), 2);
pointa = zeros(length(Theta0), 2);
pointb = zeros(length(Theta0), 2);
pointc = zeros(length(Theta0), 2);
pointd = zeros(length(Theta0), 2);
pointe = zeros(length(Theta0), 2);
pointf = zeros(length(Theta0), 2);
pointg = zeros(length(Theta0), 2);
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
% Counters
error_count = 0;
iteration_count = zeros(length(Theta0), 5);
% Calculate all angles
tic
for N = 1:length(Theta0)
    if N==1
    [Angles(N,:),error_count,iteration_count(N,:)] = JansenV2(Theta0(N),[150 210 80 180 270 180 290 250 260 300]*(pi/180), error_count);
    else
    %
    [Angles(N,:),error_count,iteration_count(N,:)] = JansenV2(Theta0(N), Angles(N-1,:),error_count); 
    end
end
toc
tic
for N = 1:length(Theta0)
    % From the angles, compute the positions
    point0(N,:) = [0, 0];
    pointc(N,:) = [-38, -7.8];
    pointa(N,:) = [a0*cos(Theta0(N)), a0*sin(Theta0(N))];
    pointb(N,:) = [a0*cos(Theta0(N))+ab*cos(Angles(N, 1)), a0*sin(Theta0(N))+ab*sin(Angles(N, 1))];
    pointd(N,:) = [-38+cd*cos(Angles(N, 4)),-7.8+cd*sin(Angles(N, 4))];
    pointe(N,:) = [-38+cd*cos(Angles(N, 4))+de*cos(Angles(N, 5)),-7.8+cd*sin(Angles(N, 4))+de*sin(Angles(N, 5))];
    pointf(N,:) = [-38+cf*cos(Angles(N, 7)),-7.8+cf*sin(Angles(N, 7))];
    pointg(N,:) = [pointe(N,1)+eg*cos(Angles(N, 10)),pointe(N,2)+eg*sin(Angles(N, 10))];
    % Open a figure to plot the linkage
    jansenfig = figure('Position',[100 100 1600 900]);
    axis([-120 120 -100 100])
    axis off
    % Creating linkage
    vertices = [point0(N,:), 0; pointa(N,:), 0; pointb(N,:), 0; pointc(N,:), 0; pointd(N,:), 0; pointe(N,:), 0; pointf(N,:), 0; pointg(N,:), 0]; 
    linkmatrix = [1 2 2 2;2 3 4 7;3 4 5 5;4 5 6 7;6 7 8 8];
    patch('Vertices',vertices, 'Faces',linkmatrix, 'FaceColor','none'); hold on;
    plot(pointg([1:N], 1), pointg([1:N],2), 'r')
    % Save frame to a movie
    M(N) = getframe(jansenfig);
    close(jansenfig)
end
toc
figure('Position',[100 100 1600 900]);
axis([-120 50 -100 100])
axis off
movie(M, 1)
if(error_count ~= 0)
    warning('Newton Raphson did not converge %d times.', error_count)
end
% 
V = VideoWriter('test', 'MPEG-4');
open(V)
writeVideo(V,M)
close(V)
y = [iteration_count([1:N],1) iteration_count([1:N],2) iteration_count([1:N],3) iteration_count([1:N],4) iteration_count([1:N],5)];
figure
bar(Theta0*(pi/180), y)
axis([ -0.5*(pi/180) 10.5*(pi/180) 0 7])
xlabel('Crank Angle (radians)')
ylabel('Number of iterations')
legend('Closure pair 1','Closure pair 2','Closure pair 4','Closure pair 3','Closure pair 5');
