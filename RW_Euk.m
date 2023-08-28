clear
clc
close all

load('x_dist_Eu_tunnel.mat')
load('r_Eu_tunnel.mat')

// %%
L = max(x_tunnel);

// %% interp
new_x = linspace(min(x_tunnel), max(x_tunnel), length(x_tunnel));   % new x values
new_r = interp1(x_tunnel,r, new_x);   % interpolate new y values  % interpolate new r values
x_tunnel = new_x;
y_tunnel = new_r;

// %% smooting
smoothing_method = 'sgolay';
// % Specify the smoothing window size or span
smoothing_window = 55;
// % Smooth the y values using the specified method and window
smoothed_y = smoothdata(y_tunnel, smoothing_method, smoothing_window);
y_tunnel = smoothed_y;
// %%
// % plot(x_tunnel, y_tunnel,'b-','linewidth',1.5)
// % hold on
// % plot(x_tunnel, -y_tunnel,'b-','linewidth',1.5)
// % hold on 
// % Define parameters
// %%
// % R = [1 0.7 0.5 0.3 0.2 0.1];
R = 1;
// % plot(x_tunnel, R*y_tunnel,'r-','linewidth',1.5)
// % hold on
// % plot(x_tunnel, -R*y_tunnel,'r-','linewidth',1.5)
// % hold on
x_init = x_tunnel(10);  % Initial x-position of the particle
K = 15000;  % Number of simulation runs
D = 0.5;
dt=0.01; %delta_t
M = length(R);
T = 1000000;
MFPT = zeros(M,1);
% Initialize variables
% h = animatedline;
FPT=zeros(K,M);
for n=1:M
    parfor k = 1:K
        if mod(k,100)==0
                k
        end
        // % Initialize variables for each simulation run
        x = x_init;  % Current x-position of the particle
        y = 0.1;  % Current y-position of the particle
        t = 0;  % Current time
        // % Define tunnel boundaries
        y_lower = @(x) -interp1(x_tunnel, R(n)*y_tunnel, x);  % Lower boundary function
        y_upper = @(x) interp1(x_tunnel, R(n)*y_tunnel, x);  % Upper boundary function
        // % Perform the random walk simulation
        while x > 0 && x < L && t<T %&& y > y_lower(x) && y < y_upper(x)
            dx = sqrt(2*D*dt)*(randn);  % Random displacement in the x-direction
            dy = sqrt(2*D*dt)*randn;  % Random displacement in the y-direction

            x = x + dx;
            y = y + dy;

            if x <= 0
                x = x - dx;
            end
            if y < y_lower(x) || y > y_upper(x)
                y = y - dy;
            end
            t = t + 1;
        end
        if x>=L
            FPT(k,n) = t;
        end
    end
end
for n=1:M
    FPT = FPT(:,n);
    FPT = FPT(FPT ~= 0);
end
histogram(FPT,20);
MFPT(n) = mean(FPT)

// % mywriter = VideoWriter('RW_E_R03_abs');
// % mywriter.FrameRate = 20;
// % open(mywriter);
// % writeVideo(mywriter,movievector);
// % close(mywriter);
// % Calculate the mean first passage time
// % meanFirstPassageTime = meanFirstPassageTime / K;

// % plot(R, MFPT,'ob-','linewidth',1.5,'MarkerFaceColor','k','markersize',5)
// % xlabel('R')
// % ylabel('Mean First passage time')
// % title('RW a Single Particle in Eukaryotic Tunnel');

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %         y_tunnel(120:140) = y_tunnel(120:140)*R(n); %% first bottleneck Euk
// %         y_tunnel(210:250) = y_tunnel(210:250)*R(n); %% zero bottleneck Euk
// %         y_tunnel(340:375) = y_tunnel(340:375)*R(n); %% second bottleneck Euk
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%