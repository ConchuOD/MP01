function [M, theta] = Jensen_Mech_caller()
% uncommented mostly - for plotting only
    OI = 15;
    IJ = 50;
    JK = 55.8;
    KL = 39.4;
    LM = 65.7;
    MN = 49;
    LN = 36.7;
    NI = 61.9;
    N1 = 39.3;
    K1 = 40.1;
    J1 = 41.5;
    
    consts = zeros(5, 2);
    theta = zeros(1,11);
    vert_leg = NaN(360,2);

    coeffs(1,:) = [J1,-IJ];
    coeffs(2,:) = [K1,-JK]; 
    coeffs(3,:) = [-N1,-NI];%NI
    coeffs(4,:) = [-KL,-LN];%KL
    coeffs(5,:) = [-LM,MN];%-LM

    theta_est = [(5*pi)/9,2*pi/3,2*pi/3,10*pi/9,11*pi/18,11*pi/18,10*pi/9,11*pi/18,pi/6,5*pi/6]';
    
    t = linspace(0, 2*pi);
    R = 1;
    circ_vec = [R*cos(t);R*sin(t)];
    
    x_rot = (OI/R)*circ_vec(1,:);
    y_rot = (OI/R)*circ_vec(2,:);
    rotator_vertices = [x_rot', y_rot'];
    rot_link = [1:length(rotator_vertices(:,1))];
    
    % Joint circles which have fixed positions
    x_0 = circ_vec(1,:);
    y_0 = circ_vec(2,:);
    
    x_1 = circ_vec(1,:)-38;
    y_1 = circ_vec(2,:)-7.8;
    
    for angle_count = 1:361
        theta(1) = (361-angle_count)*pi/180;
        
        consts(1,:) = [OI*cos(theta(1))+38, OI*sin(theta(1))+7.8];
        [theta(2), theta(3)] = Jensen_Mech_Angles(coeffs(1,:), consts(1,:), theta_est(1:2));
    
        consts(2,:) = [J1*cos(theta(2)), J1*sin(theta(2))]; 
        [theta(4), theta(5)] = Jensen_Mech_Angles(coeffs(2,:), consts(2,:), theta_est(3:4));
    
        consts(3,:) = [OI*cos(theta(1))+38, OI*sin(theta(1))+7.8];%all changed from negative
        [theta(7), theta(8)] = Jensen_Mech_Angles(coeffs(3,:), consts(3,:),theta_est(6:7));
    
        consts(4,:) = [-N1*cos(theta(7))-K1*cos(theta(4)), -N1*sin(theta(7))-K1*sin(theta(4))];%t7 pos 
        [theta(6), theta(11)] = Jensen_Mech_Angles(coeffs(4,:), consts(4,:),[theta_est(5);theta_est(10)]);
    
        consts(5,:) = [-LN*cos(theta(11)), -LN*sin(theta(11))];%pos LN
        [theta(9), theta(10)] = Jensen_Mech_Angles(coeffs(5,:), consts(5,:),theta_est(8:9));
    
        % Position vectors of joints
        % joint 0
        j_0 = [0;0];
        % joint 1
        j_1 = [-38;-7.8];
        % joint i
        j_i = [OI*cos(theta(1));OI*sin(theta(1))];
        % joint j
        j_j = [-38+J1*cos(theta(2));-7.8+J1*sin(theta(2))];
        % joint k
        j_k = [-38+K1*cos(theta(4));-7.8+K1*sin(theta(4))];
        % joint l
        j_l = [j_k(1)-KL*sin(pi/2 -theta(6));j_k(2)-KL*cos(pi/2 -theta(6))];
        % joint n
        j_n = [OI*cos(theta(1))+NI*cos(theta(8));OI*sin(theta(1))+NI*sin(theta(8))];
        % joint m
        j_m = [j_n(1)-MN*cos(theta(10));j_n(2)-MN*sin(theta(10))];
    
        % Position vectors of corners of bars ha
        % Bar 0_i
        
        
        %Joint circles which change location
        x_i = j_i(1)+circ_vec(1,:);
        y_i = j_i(2)+circ_vec(2,:);
        
        x_j = j_j(1)+circ_vec(1,:);
        y_j = j_j(2)+circ_vec(2,:);
        
        x_k = j_k(1)+circ_vec(1,:);
        y_k = j_k(2)+circ_vec(2,:);
        
        x_l = j_l(1)+circ_vec(1,:);
        y_l = j_l(2)+circ_vec(2,:);
        
        x_m = j_m(1)+circ_vec(1,:);
        y_m = j_m(2)+circ_vec(2,:);
        
        x_n = j_n(1)+circ_vec(1,:);
        y_n = j_n(2)+circ_vec(2,:);
        
        % Rotator circles
        rx_1 = 7.5*cos(theta(1)+pi/2)+circ_vec(1,:);
        ry_1 = 7.5*sin(theta(1)+pi/2)+circ_vec(2,:);
        
        rx_2 = 7.5*cos(theta(1)+3*pi/2)+circ_vec(1,:);
        ry_2 = 7.5*sin(theta(1)+3*pi/2)+circ_vec(2,:);
        
        vert_leg(angle_count,:) = j_m';
        face_leg = [1:angle_count];
    
        vertices = [j_0';j_1';j_i';j_j';j_k';j_l';j_m';j_n'];
        faces = [[2 4 5 NaN];[2 5 6 8];[6 7 8 NaN];[2 8 3 4];[1 3 NaN NaN]];
        
        triangle_faces = [[2 4 5 NaN];[6 7 8 NaN]];
        
        handlefig = figure('Position', [100 100 850 600]); axis square;xlim([-120 40]);ylim([-120 40]);
        
        patch('Vertices', rotator_vertices, 'Faces', rot_link, 'FaceColor', [0.85,0.85,0.85]);
        patch(rx_1, ry_1, 'w');
        patch(rx_2, ry_2, 'w');
        patch('Vertices', vertices, 'Faces', faces, 'FaceColor', [1,1,1], 'EdgeColor', [0,0,0], 'FaceAlpha', 0);
        patch('Vertices', vertices, 'Faces', triangle_faces, 'FaceColor', [0.5,0.5,0.5], 'EdgeColor', [0,0,0]);
        patch('Vertices', vert_leg, 'Faces', face_leg, 'FaceColor', [1,1,1], 'EdgeColor', [1,0,0], 'FaceAlpha', 0);
        patch([2.1213*cos(3*pi/4+theta(1)),2.1213*cos(5*pi/4+theta(1)), j_i(1)+2.1213*cos(7/4*pi+theta(1)), j_i(1)+2.1213*cos(theta(1)+pi/4)], [2.1213*sin(3*pi/4+theta(1)),2.1213*sin(5*pi/4+theta(1)), j_i(2)+2.1213*sin(7/4*pi+theta(1)), j_i(2)+2.1213*sin(theta(1)+pi/4)], 'green');
        patch(x_0, y_0, 'g');
        patch([j_j(1)+2.1213*cos(pi/4+theta(3)),j_j(1)+2.1213*cos(theta(3)-pi/4), j_i(1)+2.1213*cos(theta(3)+5*pi/4),j_i(1)+2.1213*cos(3*pi/4+theta(3))], [j_j(2)+2.1213*sin(pi/4+theta(3)),j_j(2)+2.1213*sin(theta(3)-pi/4),j_i(2)+2.1213*sin(theta(3)+5*pi/4), j_i(2)+2.1213*sin(3*pi/4+theta(3))], 'green');
        patch([j_n(1)+2.1213*cos(theta(8)+pi/4), j_n(1)+2.1213*cos(theta(8)-pi/4), j_i(1)+2.1213*cos(theta(8)-3*pi/4), j_i(1)+2.1213*cos(theta(8)-5*pi/4)], [j_n(2)+2.1213*sin(theta(8)+pi/4), j_n(2)+2.1213*sin(theta(8)-pi/4), j_i(2)+2.1213*sin(theta(8)-3*pi/4), j_i(2)+2.1213*sin(theta(8)-5*pi/4)], 'green');
        patch([j_k(1)+2.1213*cos(theta(6)+pi/4), j_k(1)+2.1213*cos(theta(6)-pi/4), j_l(1)+2.1213*cos(theta(6)+5*pi/4), j_l(1)+2.1213*cos(theta(6)+3*pi/4)], [j_k(2)+2.1213*sin(theta(6)+pi/4), j_k(2)+2.1213*sin(theta(6)-pi/4), j_l(2)+2.1213*sin(theta(6)+5*pi/4), j_l(2)+2.1213*sin(theta(6)+3*pi/4)], 'green');
        patch([j_1(1)+2.1213*cos(theta(7)+pi/4), j_1(1)+2.1213*cos(theta(7)-pi/4), j_n(1)+2.1213*cos(theta(7)+5*pi/4), j_n(1)+2.1213*cos(theta(7)+3*pi/4)], [j_1(2)+2.1213*sin(theta(7)+pi/4), j_1(2)+2.1213*sin(theta(7)-pi/4), j_n(2)+2.1213*sin(theta(7)+5*pi/4), j_n(2)+2.1213*sin(theta(7)+3*pi/4)], 'green');

        patch(x_1, y_1, 'g');
        patch(x_i, y_i, 'g');
        patch(x_j, y_j, 'g');
        patch(x_k, y_k, 'g');
        patch(x_l, y_l, 'g');
        patch(x_m, y_m, 'g');
        patch(x_n, y_n, 'g');
        
        M(angle_count) = getframe(handlefig);
        close(handlefig)
        
        % Assign the current calculated approximations of the angles to the next estimate, assuming
        % that the step size is relatively small, this will decreasse the
        % number of iterations of the Newton-Raphson method required for
        % the solution to converge.
        theta_est = theta(2:11)';
    end
end