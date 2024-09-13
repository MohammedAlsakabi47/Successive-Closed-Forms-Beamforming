close all
clear all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------------- Inputs ------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = 16; %------------------------------------------------Number of elements
theta = [-60, -20, 20, 45, 60, 80]; %------------------------------------Directions
freq = 1*10^9; %---------------------------------------Operating frequency
lambda = (3*10^8)/freq; %--------------------------------------------Lambda
d=lambda/2;
WF = sqrt(diag([1 1 1 1 1])); %---------------------------Weighting factors
beamDirac = -90:90; %------------------------------------------Beam range
threshold = 5*10^-16;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------ Steering vectors -------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


steering = 2*pi*d/lambda; %------------------------------------Phase Factor
for L = 1 : length(theta)

    for k = 1 : M
        %V(L,k)= (sind(theta(L)))^2 * exp((-1)*steering*sind(theta(L))*(k-1)*i); %---------Transmitting steering vector: tan square * delay
        V(L,k)= exp((-1)*steering*sind(theta(L))*(k-1)*1i); %---------Transmitting steering vector: tan square * delay
    end

end

V = transpose(V);
S = (V)*(ctranspose(V));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------- EC -------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
[eigVector,eigValue] = eig(S);
eigValue = diag(eigValue);
Index = find(eigValue == min(eigValue));
x_EC = eigVector(:,Index);

Null_Power = (ctranspose(x_EC))*S*x_EC;

wheighted_Nulls_EC = Null_Power*WF;
wheighted_Nulls_EC = diag(wheighted_Nulls_EC);

elapsedTime_EC = toc



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------- CMC -------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%In CMC the signals x_4, x_8 and x_16 are fixed for testing purpose 
%Each M should have identical threshold for different # of directions



R = [real(S) -imag(S); imag(S) real(S)]; %----------Step 1 of the algorithm
lambda_fac = 0.5;
R_dash = 2*(R + lambda_fac*eye(length(R)));


x=randn(M,1)+1i*randn(M,1); %---------------------------Initialization of x

% %%%%% fixed x %%%%%
% x_16 = transpose([0.1732 + 0.4980i -0.5055 + 2.7891i -1.1933 + 0.7276i 0.6470 - 0.7731i -0.3536 + 0.8366i 0.0464 - 1.1283i -0.7929 - 1.4245i -1.5505 + 0.7174i 0.1716 - 0.7779i -0.0621 + 0.3160i 1.1990 + 1.4065i 0.8017 + 0.4011i 1.0533 + 0.9297i -0.7489 - 1.6058i -0.9363 + 0.6615i -1.2691 + 2.1385i]);
% 
% 
% 
% x_4 =[  1.1921 + 1.0205i 
%        -1.6118 + 0.8617i
%        -0.0245 + 0.0012i
%        -1.9488 - 0.0708i];
% 
% 
% x_8 =[ -0.2938 - 0.1765i
%        -0.8479 + 0.7914i
%        -1.1201 - 1.3320i
%         2.5260 - 2.3299i
%         1.6555 - 1.4491i
%         0.3075 + 0.3335i
%        -1.2571 + 0.3914i
%        -0.8655 + 0.4517i];
% 
% 
%   
% 
% if M==16
%     x = x_16;
% else if M == 8
%         x = x_8;
%     else
%         x = x_4; 
%     end
% end

        


B=zeros(M,2*M); %---------------------------------------Setting of B matrix
frst_x = x;
x_old = x;
x_new = x_old;

x=randn(M,1)+1i*randn(M,1);
xrand = x;


tic
n = 1;
counter = n;
while n >= 1
    B=zeros(M,2*M);
    for ii=1:M
        for jj=1:M
            if ii==jj
                B(ii,jj)=cos(angle(x(ii)));
                B(ii,jj+M)=sin(angle(x(ii)));
            end
        end
    end
    
    s = (inv(R_dash))*(transpose(B))*(inv(B*(inv(R_dash))*(transpose(B))))*(ones(length(R_dash)/2,1)); %Step 3
    
    for index = 1:M
        s_complex = s(index)+1i*s(index+M);
        x_new(index,:) = exp(1i*(angle(s_complex)));
    end

    x=x_new;
    fx(n,:) = ctranspose(x_new)*S*x_new;
    fd(n) =(abs(ctranspose(x_new)*S*x_new));

    if n>1
        if abs(abs(fd(n))-abs(fd(n-1)))
            fd_final = fd(n);
            fd_count(n) = fd_final;
        end
    else
        fd_final = fd(n);
        fd_count(n) = fd_final;
    end
    
    if fd_final <= threshold
        break
    end
    
    n = n+1;
    counter = counter+1;

end

elapsedTime_CMC = toc
fd_final=min(fd);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------ plot CMC ---------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v=zeros(M,1);
for ii=1:length(beamDirac)
    for k=1:M
        v(k)= exp((-1i)*steering*sind(beamDirac(ii))*(k-1));
    end
    F1(ii)=abs(v'*x).^2;
end

figure,CMC_plot = plot(beamDirac,10*log10(F1/max(F1)), 'LineWidth',3);
hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------ plot EC ----------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v=zeros(M,1);
for ii=1:length(beamDirac)
    for k=1:M
        v(k)= exp((-1i)*steering*sind(beamDirac(ii))*(k-1));
    end
    F2(ii)=abs(v'*x_EC).^2; 
end

EC_plot = plot(beamDirac,10*log10(F2/max(F2)));
xlim([beamDirac(1) beamDirac(length(beamDirac))])


EC_result = mag2db(abs(Null_Power))
CMC_result = mag2db(fd_final)


for k = 1 : length(theta)
    
    plot([theta(k) theta(k)], [-350 0],'k--')
    numstg = num2str(theta(k));
    txt = [num2str(theta(k)),char(176), '\rightarrow'];
    text(theta(k),-20,txt,'HorizontalAlignment','right','FontSize',12)

end

legend([CMC_plot EC_plot], {'CMC','EC'},'Location','east')
 xlabel('Angle (degree)')
 ylabel('Array Factor (dB)')
grid on
