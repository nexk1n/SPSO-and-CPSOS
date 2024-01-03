% Implement Particle Swarm Optimization (PSO) algorithm
clear;
close all;

iteration_best_fitness_SPSO = zeros(1000, 50);
iteration_best_fitness_CPSOS = zeros(1000, 50);

for n = 1:50
    % Use Rosenbrock as the objective function
    % Random initialization of 100 particles in the range [-10, 10] dimension 10
    partcles_SPSO = rand(100, 10) * 20 - 10;
    partcles_CPSOS = rand(100, 10) * 20 - 10;
    % Chaotic initialization of 100 particles
    ch = rand(100, 10);

    for i = 1:100
        partcles_CPSOS(i, :) = sin(pi * ch(i, :));
    end

    % Initialize velocity of particles
    velocity_SPSO = zeros(100, 10);
    velocity_CPSOS = zeros(100, 10);

    % Calculate the fitness of each particle
    fitness_SPSO = zeros(100, 1);
    fitness_CPSOS = zeros(100, 1);

    for i = 1:9
        fitness_SPSO = fitness_SPSO + 100 * (partcles_SPSO(:, i + 1) - partcles_SPSO(:, i) .^ 2) .^ 2 + (partcles_SPSO(:, i) - 1) .^ 2;
        fitness_CPSOS = fitness_CPSOS + 100 * (partcles_CPSOS(:, i + 1) - partcles_CPSOS(:, i) .^ 2) .^ 2 + (partcles_CPSOS(:, i) - 1) .^ 2;
    end

    % Initialize the best position of each particle
    best_position_SPSO = partcles_SPSO;
    best_position_CPSOS = partcles_CPSOS;
    % Initialize the best fitness of each particle
    best_fitness_SPSO = fitness_SPSO;
    best_fitness_CPSOS = fitness_CPSOS;
    % Initialize the best position of the whole swarm
    [global_best_fitness_SPSO, global_best_index_SPSO] = min(fitness_SPSO);
    global_best_position_SPSO = partcles_SPSO(global_best_index_SPSO, :);
    [global_best_fitness_CPSOS, global_best_index_CPSOS] = min(fitness_CPSOS);
    global_best_position_CPSOS = partcles_CPSOS(global_best_index_CPSOS, :);
    % Initialize the parameters
    c1_SPSO = 2;
    c2_SPSO = 2;
    w_SPSO = 0.9;
    c1_CPSOS = 1;
    c2_CPSOS = 1;
    w_CPSOS = 0.9;

    % Record the best fitness and position of all iterations
    optimal_fitness_SPSO = global_best_fitness_SPSO;
    optimal_position_SPSO = global_best_position_SPSO;
    optimal_fitness_CPSOS = global_best_fitness_CPSOS;
    optimal_position_CPSOS = global_best_position_CPSOS;

    iteration_best_fitness_SPSO(1, n) = global_best_fitness_SPSO;
    iteration_best_fitness_CPSOS(1, n) = global_best_fitness_CPSOS;

    for i = 1:999
        % Update the velocity of each particle
        velocity_SPSO = w_SPSO * velocity_SPSO + c1_SPSO * rand(100, 10) .* (best_position_SPSO - partcles_SPSO) + c2_SPSO * rand(100, 10) .* (global_best_position_SPSO - partcles_SPSO);
        velocity_CPSOS = w_CPSOS * velocity_CPSOS + c1_CPSOS * rand(100, 10) .* (best_position_CPSOS - partcles_CPSOS) + c2_CPSOS * rand(100, 10) .* (global_best_position_CPSOS - partcles_CPSOS);
        % Update the position of each particle
        partcles_SPSO = partcles_SPSO + velocity_SPSO;
        % Calculate relative fitness of each particle
        lambda_CPSOS = fitness_CPSOS - sum(fitness_CPSOS) ./ 100;
        lambda_CPSOS = exp(lambda_CPSOS);

        for j = 1:100

            if lambda_CPSOS(j) < rand(1)
                r = rand(1) * 2 - 1;
                partcles_CPSOS(j, :) = abs(partcles_CPSOS(j, :) - velocity_CPSOS(j, :)) * exp(r) * cos(2 * pi * r) + global_best_position_CPSOS;
            else
                partcles_CPSOS(j, :) = partcles_CPSOS(j, :) + velocity_CPSOS(j, :);
            end

        end

        % Update the fitness of each particle
        fitness_SPSO = zeros(100, 1);
        fitness_CPSOS = zeros(100, 1);

        for j = 1:9
            fitness_SPSO = fitness_SPSO + 100 * (partcles_SPSO(:, j + 1) - partcles_SPSO(:, j) .^ 2) .^ 2 + (partcles_SPSO(:, j) - 1) .^ 2;
            fitness_CPSOS = fitness_CPSOS + 100 * (partcles_CPSOS(:, j + 1) - partcles_CPSOS(:, j) .^ 2) .^ 2 + (partcles_CPSOS(:, j) - 1) .^ 2;
        end

        % Update the best position of each particle
        best_position_SPSO(fitness_SPSO < best_fitness_SPSO, :) = partcles_SPSO(fitness_SPSO < best_fitness_SPSO, :);
        best_fitness_SPSO(fitness_SPSO < best_fitness_SPSO) = fitness_SPSO(fitness_SPSO < best_fitness_SPSO);
        best_position_CPSOS(fitness_CPSOS < best_fitness_CPSOS, :) = partcles_CPSOS(fitness_CPSOS < best_fitness_CPSOS, :);
        best_fitness_CPSOS(fitness_CPSOS < best_fitness_CPSOS) = fitness_CPSOS(fitness_CPSOS < best_fitness_CPSOS);
        % Update the best position of the whole swarm
        [global_best_fitness_SPSO, global_best_index_SPSO] = min(best_fitness_SPSO);
        global_best_position_SPSO = best_position_SPSO(global_best_index_SPSO, :);
        [global_best_fitness_CPSOS, global_best_index_CPSOS] = min(best_fitness_CPSOS);
        global_best_position_CPSOS = best_position_CPSOS(global_best_index_CPSOS, :);
        % Update the parameters
        w_SPSO = w_SPSO - 0.5/1000;
        w_CPSOS = w_CPSOS - 0.5/1000;
        c1_CPSOS = 0.5 + 0.5 * exp(-i / 500) + 1.4 * sin(i) / 30;
        c2_CPSOS = 1 + 1.4 * (1 - exp(-i / 500)) + 1.4 * sin(i) / 30;
        % Record the best fitness and position of all iterations
        if global_best_fitness_SPSO < optimal_fitness_SPSO
            optimal_fitness_SPSO = global_best_fitness_SPSO;
            optimal_position_SPSO = global_best_position_SPSO;
        end

        if global_best_fitness_CPSOS < optimal_fitness_CPSOS
            optimal_fitness_CPSOS = global_best_fitness_CPSOS;
            optimal_position_CPSOS = global_best_position_CPSOS;
        end

        iteration_best_fitness_SPSO(i + 1, n) = global_best_fitness_SPSO;
        iteration_best_fitness_CPSOS(i + 1, n) = global_best_fitness_CPSOS;
    end

    % Print the optimal fitness and position of SPSO and CPSOS
    fprintf('The optimal fitness of SPSO is %f\n', optimal_fitness_SPSO);
    fprintf('The optimal position of SPSO is\n');
    disp(optimal_position_SPSO);
    fprintf('The optimal fitness of CPSOS is %f\n', optimal_fitness_CPSOS);
    fprintf('The optimal position of CPSOS is\n');
    disp(optimal_position_CPSOS);

end

% Plot the average best fitness of all iterations of SPSO and CPSOS in one figure
figure(1);
plot(mean(iteration_best_fitness_SPSO, 2), 'r');
hold on;
plot(mean(iteration_best_fitness_CPSOS, 2), 'b');
xlabel('Iteration');
ylabel('Average best fitness');
legend('SPSO', 'CPSOS');
title('Average best fitness of all iterations of SPSO and CPSOS (Rosenbrock)');
% print average final best fitness
fprintf('The average final best fitness of SPSO is %f\n', mean(iteration_best_fitness_SPSO(1000, :)));
fprintf('The average final best fitness of CPSOS is %f\n', mean(iteration_best_fitness_CPSOS(1000, :)));
% print standard deviation of final best fitness
fprintf('The standard deviation of final best fitness of SPSO is %f\n', std(iteration_best_fitness_SPSO(1000, :)));
fprintf('The standard deviation of final best fitness of CPSOS is %f\n', std(iteration_best_fitness_CPSOS(1000, :)));
