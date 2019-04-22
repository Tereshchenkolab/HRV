function [resp_points_US,resp_points_DS,resp_points_RA] = resp_points(xi,q_points_x,r_points_x,s_points_x,fs,samp8)


I_US=[];
I_DS=[];
phi_R=[];
resp_points_US=[];
resp_points_DS=[];
resp_points_RA=[];

for i=1:length(r_points_x)

    n_qr=xi(q_points_x(i):r_points_x(i));
    n_rs=xi(r_points_x(i):s_points_x(i));

    [n_U_val,n_U_ind] = max(abs(diff(n_qr)));
    [n_D_val,n_D_ind] = max(abs(diff(n_rs)));


    %Up
    U_y=n_qr(max(1,n_U_ind-round(samp8/2)):min(n_U_ind+round(samp8/2),length(n_qr)));
    U_x=1:length(U_y);

    % Formulating in matrix for solving for least squares estimate
    Y1 = U_y.';
    X1 = [U_x.' ones(1,length(U_x)).'];
%     alpha = inv(X'*X)*X'*Y; % solving for m and c
    alpha1 = inv(X1'*X1)*X1'*Y1'; % solving for m and c
    % constructing the straight line using the estimated slope and constant
    yEst1 = alpha1(1)*U_x + alpha1(2);

%     figure
%     plot(U_x,U_y,'r.')
%     hold on
%     plot(U_x,yEst,'b')
%     legend('observations', 'estimated straight line')
%     grid on
%     ylabel('observations')
%     xlabel('x axis')
%     title('least squares straight line fit')

%     figure;
%     plot(n_qr);hold on;plot(U_x+max(1,n_U_ind-round(samp8/2)-1),U_y,'ro');plot(U_x+max(1,n_U_ind-round(samp8/2)-1),yEst1,'k');hold off


    %Down
    D_y=n_rs(max(1,n_D_ind-round(samp8/2)):min(n_D_ind+round(samp8/2),length(n_rs)));
    D_x=1:length(D_y);

    % Formulating in matrix for solving for least squares estimate
    Y2 = D_y.';
    X2 = [D_x.' ones(1,length(D_x)).'];
%     alpha = inv(X'*X)*X'*Y; % solving for m and c
    alpha2 = inv(X2'*X2)*X2'*Y2'; % solving for m and c
    % constructing the straight line using the estimated slope and constant
    yEst2 = alpha2(1)*D_x + alpha2(2);

%     figure
%     plot(U_x,U_y,'r.')
%     hold on
%     plot(U_x,yEst,'b')
%     legend('observations', 'estimated straight line')
%     grid on
%     ylabel('observations')
%     xlabel('x axis')
%     title('least squares straight line fit')

%     figure;
%     plot(n_rs);hold on;plot(D_x+max(1,n_D_ind-round(samp8/2)-1),D_y,'ro');plot(D_x+max(1,n_D_ind-round(samp8/2)-1),yEst2,'k');hold off


    %Slopes of lines
    I_US(i)=alpha1(1);
    I_DS(i)=alpha2(1);

    %R-wave angle
    phi_R(i)=atan(abs((I_US(i)-I_DS(i))/0.4*(6.25+(I_US(i)*I_DS(i))))); 


    %Respiration Upper
    resp_points_US(i,1)=r_points_x(i);
    resp_points_US(i,2)=I_US(i);

    %Respiration Lower
    resp_points_DS(i,1)=r_points_x(i);
    resp_points_DS(i,2)=I_DS(i);

    %Respiration R_angle
    resp_points_RA(i,1)=r_points_x(i);
    resp_points_RA(i,2)=phi_R(i);

end

end