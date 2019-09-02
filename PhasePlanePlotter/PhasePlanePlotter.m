clc
clear
close all

disp('Phase Plane Plotter Software');
disp('Created by Masoud Shakoorian Fard');
disp(' ');

% Define x1 and x2 Boundaries
disp('---- Boundary Definition Part ----');
disp('Please Enter x1 and x2 Boundary in the form [x1_min x1_max x2_min x2_max]');
disp('The Default Value is [-0.5 0.5 -0.5 0.5]. To choose Default Value Enter 0');
x_bound = input('Boundaries: ');

if x_bound == 0
    x_bound = [-0.5 0.5 -0.5 0.5];
end
[x1, x2] = meshgrid(linspace(x_bound(1), x_bound(2), 20), linspace(x_bound(3), x_bound(4), 20));

disp(' ');

% Define State Equations
disp('---- State Equation Definition Part ----');
eqns_num = 2;
while 1
%     eqns_num = input('Please Enter Number of State Equations: ');

    if eqns_num == 1
        disp('Please Enter State Equation in this form [Example: -x^3]:');
        xdot_str = input('xdot = ', 's');
        xdot = eval(vectorize(xdot_str));
        disp(' ');
        break;

    elseif  eqns_num == 2
        for i=1:eqns_num
            if i == 1
                disp('Please Enter First State Equation in this form [Example: -x1 - 2*x2*x1^2 + x2]:');
                x1dot_str = input('x1dot = ', 's');
                x1dot = eval(vectorize(x1dot_str));
                xdot1 = replace(x1dot_str,["x1" , "x2"],["x(1)","x(2)"]);
                disp(' ');

            elseif i == 2
                disp('Please Enter Second State Equation in this form [Example: -x1 - x2]:');
                x2dot_str = input('x2dot = ', 's');
                x2dot = eval(vectorize(x2dot_str));
                xdot2 = replace(x2dot_str,["x1" , "x2"],["x(1)","x(2)"]);
                disp(' ');
            end
        end
        break;

    else
        disp('Not Supported. Please Try Again!');
        disp(' ');
    end
end

state_function = "@(t,x) [" + xdot1 + ";" + xdot2 + "]";
state_function = str2func(state_function);

% To Find Eq. Points
state_function_e = "@(x) [" + xdot1 + ";" + xdot2 + "]";
state_function_e = str2func(state_function_e);

% Vector Field Plot
linewidth = 1; %Linewidth for Trajectories and Lyapanov Functions

figure('Name','Phase Plane Plotter','NumberTitle','off');
quiver(x1,x2,x1dot,x2dot);
hold on;
grid on;
% axis equal;
xlim([x_bound(1) x_bound(2)]);
ylim([x_bound(3) x_bound(4)]);
title('Phase Plane Plot');
xlabel('x_1');
ylabel('x_2');
box on
set(gca,'fontsize',16,'FontWeight','bold')
figureHandle = gcf;
set(findall(figureHandle,'type','text'),'FontName','Times New Roman','fontSize',16,'fontWeight','bold','FontAngle','oblique')
set(gca,'FontName','Times New Roman')
set(gca,'fontSize',20)

% Eq. Points Plot
disp('Eq. Points: ');
x10_e = linspace(x_bound(1), x_bound(2), 20);
x20_e = linspace(x_bound(3), x_bound(4), 20);

eq_temp = [1000000 1000000];
k = 1;
for i=1:20
    options = optimoptions('fsolve','Display','none');
    [eq, ~, exit_flag, ~] = fsolve(state_function_e,[x10_e(i) x20_e(i)], options);
    if exit_flag > 0
        plot(eq(1), eq(2), 'bo', 'Linewidth', 4*linewidth);
        hold on;

        if norm(eq - eq_temp , 2) > 0.01
            eq_points(k,:) = eq;
            k = k + 1;
            disp(eq);
        end
        eq_temp = eq;
    end
end
disp(' ');

% Trajectory and Lyapanov Function Plot
disp('---- Trajectory and Lyapanov Function Plot ----');
disp('Please Choose a Condition: [Only Trajectory(1), Trajectory and Lyapanov Function(2), Only Lyapanov Function(3), Phase Portrait(4), Phase Portrait and Lyapanov Function(5)]')
disp('The Default Condition is (1)')
condition = input('Condition: ');
disp(' ')
if ((condition ~= 1) && (condition ~= 2) && (condition ~= 3) && (condition ~= 4) && (condition ~= 5))
    condition = 1;
end

if condition == 1
    disp('Please Choose Trajectory Plot Method: [Choose Initial Point(s) From Plot(1), Enter Initial Point(s)(2)]')
    disp('The Default Trajectory Plot Method is (1)')
    plot_method = input('Method: ');
    disp(' ')
    if ((plot_method ~= 1) && (plot_method ~= 2))
        plot_method = 1;
    end
    
    while 1
        if plot_method == 1
            disp('Please Choose Initial Point(s) from Figure and then Press Enter')
            disp(' ')
            clear x10 x20;
            [x10, x20] = getpts;
            for i=1:length(x10)
                clear x t
                [t,x] = ode45(state_function,[0 20],[x10(i) x20(i)]);
                plot(x(:,1),x(:,2),'k', 'Linewidth', linewidth);
                hold on
            end
        
        elseif plot_method == 2
            disp('Please Enter Initial Point(s) [Example: [0 1 2 3]]: ')
            disp(' ')
            clear x10 x20;
            x10 = input('x1 Initial Values: ');
            x20 = input('x2 Initial Values: ');
            disp(' ')
            for i=1:length(x10)
                clear x t
                [t,x] = ode45(state_function,[0 20],[x10(i); x20(i)]);
                plot(x(:,1),x(:,2),'k', 'Linewidth', linewidth);
                hold on                
            end            
        end  
    end
    
elseif condition == 2
    while 1
        traj_lyap = input('[Trajectory - Choose Points(1)] or [Trajectory - Enter Points(2)] or [Lyapanov Function(3)]: [Default: 1] ');
        disp(' ');
        if ((traj_lyap ~= 1) && (traj_lyap ~= 2) && (traj_lyap ~= 3))
            traj_lyap = 1;
        end
        
        if traj_lyap == 1
            disp('Please Choose Initial Point(s) from Figure and then Press Enter')
            disp(' ')
            clear x10 x20;
            [x10, x20] = getpts;
            for i=1:length(x10)
                clear x t
                [t,x] = ode45(state_function,[0 20],[x10(i) x20(i)]);
                plot(x(:,1),x(:,2),'k', 'Linewidth', linewidth);
                hold on
            end
            
        elseif traj_lyap == 2
            disp('Please Enter Initial Point(s) [Example: [0 1 2 3]]: ')
            disp(' ')
            clear x10 x20;
            x10 = input('x1 Initial Values: ');
            x20 = input('x2 Initial Values: ');
            disp(' ')
            for i=1:length(x10)
                clear x t
                [t,x] = ode45(state_function,[0 20],[x10(i); x20(i)]);
                plot(x(:,1),x(:,2),'k', 'Linewidth', linewidth);
                hold on                
            end   
            
        elseif traj_lyap == 3
            disp("Please Enter Lyapanov Function [Example: x^2 + y^2 - 0.25]")
            lyapanov_function = vectorize(input('Lyapanov Function: ' , 's'));
            lyapanov_function = "@(x,y) " + lyapanov_function;
            lyapanov_function = str2func(lyapanov_function);
            fimplicit(lyapanov_function, 'r', 'Linewidth', linewidth);
            hold on
        end
        disp(' ');
    end
    
elseif condition == 3
    disp("Please Enter Lyapanov Function [Example: x^2 + y^2 - 0.25]")
    while 1
        lyapanov_function = vectorize(input('Lyapanov Function: ' , 's'));
        lyapanov_function = "@(x,y) " + lyapanov_function;
        lyapanov_function = str2func(lyapanov_function);
        fimplicit(lyapanov_function, 'r', 'Linewidth', linewidth);
        hold on
    end
    
elseif condition == 4
    x10_p = linspace(x_bound(1), x_bound(2), 10);
    x20_p = linspace(x_bound(3), x_bound(4), 10);
    disp("Please Wait...")
    for i=1:length(x10_p)
        for j=1:length(x20_p)
            clear x t
            [t,x] = ode45(state_function,[0 20],[x10_p(i); x20_p(j)]);
            plot(x(:,1),x(:,2),'k', 'Linewidth', linewidth);
            hold on
        end
    end
    
    for i=1:k-1
        plot(eq_points(i,1), eq_points(i,2), 'bo', 'Linewidth', 4*linewidth);
        hold on;
    end
    disp("Finished!")
    
elseif condition == 5
    x10_p = linspace(x_bound(1), x_bound(2), 10);
    x20_p = linspace(x_bound(3), x_bound(4), 10);
    disp("Please Wait...")
    for i=1:length(x10_p)
        for j=1:length(x20_p)
            clear x t
            [t,x] = ode45(state_function,[0 20],[x10_p(i); x20_p(j)]);
            plot(x(:,1),x(:,2),'k', 'Linewidth', linewidth);
            hold on
        end
    end
    
    for i=1:k-1
        plot(eq_points(i,1), eq_points(i,2), 'bo', 'Linewidth', 4*linewidth);
        hold on;
    end
    disp("Finished!")
    
    disp("Please Enter Lyapanov Function [Example: x^2 + y^2 - 0.25]")
    while 1
        lyapanov_function = vectorize(input('Lyapanov Function: ' , 's'));
        lyapanov_function = "@(x,y) " + lyapanov_function;
        lyapanov_function = str2func(lyapanov_function);
        fimplicit(lyapanov_function, 'r', 'Linewidth', linewidth);
        hold on
    end    
end