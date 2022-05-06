function varargout = GA_TSP(varargin)

t=[40 50;160 258;100 450;200 550;300 650;350 550;290 310;370 260;470 270;500 513;510 170;140 70;400 100;267 225;210 390;380 420;50 260;280 530;200 130;100 600;260 410;540 380;70 580;120 360;300 20]; %24个点,第25个点事origin
save t.mat t
load t.mat

narginchk(0, 9);
num_custom = 24; num_dims = 2;
custom = 10*rand(num_custom, num_dims);
pop_size = 100; num_iter = 500; mutate_rate = 0.8;
show_progress = 1; show_results = 0;

% Process Inputs
custom_flag = 0; option_flag = 0;
for var = varargin
    if option_flag
        if ~isfloat(var{1}), error(['Invalid value for option ' upper(option)]); end
        switch option
            case 'popsize', pop_size = 4*ceil(real(var{1}(1))/4); option_flag = 0;
            case 'mrate', mutate_rate = min(abs(real(var{1}(1))), 1); option_flag = 0;
            case 'numiter', num_iter = round(real(var{1}(1))); option_flag = 0;
            otherwise, error(['Invalid option ' upper(option)])
        end
    elseif ischar(var{1})
        switch lower(var{1})
            case '-noplot', show_progress = 0;
            case '-results', show_results = 1;
            otherwise, option = lower(var{1}); option_flag = 1;
        end
    elseif isfloat(var{1})
        if custom_flag, error('custom or NUM_custom may be specified, but not both'); end
        if length(var{1}) == 1
            num_custom = round(real(var{1}));
            if num_custom < 2, error('NUM_custom must be an integer greater than 1'); end
            custom = 10*rand(num_custom, num_dims); custom_flag = 1;
        else
            custom = real(var{1});
            [num_custom, num_dims] = size(custom); custom_flag = 1;
            if or(num_custom < 2, num_dims ~= 2)
                error('custom must be an Nx2 matrix of floats, with N > 1')
            end
        end
    else
        error('Invalid input argument.')
    end
end

% Construct the Distance Matrix
mat3d1 = reshape(custom, 1, num_custom, num_dims);
mat3d2 = reshape(custom, num_custom, 1, num_dims);
dist_matx = sqrt(sum((mat3d1(ones(num_custom, 1), :, :) - mat3d2(:, ones(num_custom, 1), :)).^2, 3));

% Plot custom and Distance Matrix in a Figure
if show_progress
    pfig = figure;
    subplot(2, 2, 1)
    plot(custom(:, 1), custom(:, 2), 'b.')
    if num_custom < 75
        for c = 1:num_custom
            text(custom(c, 1), custom(c, 2), [' ' num2str(c)], 'Color', 'k', 'FontWeight', 'b')
        end
    end
    title([num2str(num_custom) ' custom'])
    subplot(2, 2, 2)
    imagesc(dist_matx)
    title('Distance Matrix')
    colormap(flipud(gray))
end

% Initialize Population
pop = zeros(pop_size, num_custom);
pop(1, :) = (1:num_custom);
for k = 2:pop_size
    pop(k, :) = randperm(num_custom);
end

display_rate = 20;
if num_iter < 50, display_rate = 2; end
fitness = zeros(1, pop_size);
best_fitness = zeros(1, num_iter);
for iter = 1:num_iter
    for p = 1:pop_size
        d = dist_matx(pop(p, 1), pop(p, num_custom));
        for city = 2:num_custom
            d = d + dist_matx(pop(p, city-1), pop(p, city));
        end
        fitness(p) = d;
    end
    [best_fitness(iter) index] = min(fitness);
    best_route = pop(index, :);

    % Plots
    if ~mod(iter, display_rate) && show_progress
        figure(pfig)
        subplot(2, 2, 3)
        route = custom([best_route best_route(1)], :);
        plot(route(:, 1), route(:, 2)', 'b.-')
        title(['Best GA Route (dist = ' num2str(best_fitness(iter)) ')'])
        subplot(2, 2, 4)
        plot(best_fitness(1:iter), 'r', 'LineWidth', 2)
        axis([1 max(2, iter) 0 max(best_fitness)*1.1])
    end

    % Genetic Algorithm Search
    pop = genetic_algorithm(pop, fitness, mutate_rate);
end

if show_progress
    figure(pfig)
    subplot(2, 2, 3)
    route = custom([best_route best_route(1)], :);
    plot(route(:, 1), route(:, 2)', 'b.-')
    title(['Best GA Route (dist = ' num2str(best_fitness(iter)) ')'])
    subplot(2, 2, 4)
    plot(best_fitness(1:iter), 'r', 'LineWidth', 2)
    title('Best Fitness')
    xlabel('Generation')
    ylabel('Distance')
    axis([1 max(2, iter) 0 max(best_fitness)*1.1])
end

if show_results
    figure(2)
    imagesc(dist_matx)
    title('Distance Matrix')
    colormap(flipud(gray))
    figure(3)
    plot(best_fitness(1:iter), 'r', 'LineWidth', 2)
    title('Best Fitness')
    xlabel('Generation')
    ylabel('Distance')
    axis([1 max(2, iter) 0 max(best_fitness)*1.1])
    figure(4)
    route = custom([best_route best_route(1)], :);
    plot(route(:, 1), route(:, 2)', 'b.-')
    for c = 1:num_custom
        text(custom(c, 1), custom(c, 2), [' ' num2str(c)], 'Color', 'k', 'FontWeight', 'b')
    end
    title(['Best GA Route (dist = ' num2str(best_fitness(iter)) ')'])
end

[not_used indx] = min(best_route);
best_ga_route = [best_route(indx:num_custom) best_route(1:indx-1)];
if best_ga_route(2) > best_ga_route(num_custom)
    best_ga_route(2:num_custom) = fliplr(best_ga_route(2:num_custom));
end
varargout{1} = custom(best_ga_route, :);
varargout{2} = best_ga_route;
varargout{3} = best_fitness(iter);

% --- subfunction: genetic algorithm
function new_pop = genetic_algorithm(pop, fitness, mutate_rate)
[p, n] = size(pop);

% Tournament Selection - Round One
new_pop = zeros(p, n);
ts_r1 = randperm(p);
winners_r1 = zeros(p/2, n);
tmp_fitness = zeros(1, p/2);
for k = 2:2:p
    if fitness(ts_r1(k-1)) > fitness(ts_r1(k))
        winners_r1(k/2, :) = pop(ts_r1(k), :);
        tmp_fitness(k/2) = fitness(ts_r1(k));
    else
        winners_r1(k/2, :) = pop(ts_r1(k-1), :);
        tmp_fitness(k/2) = fitness(ts_r1(k-1));
    end
end

% Tournament Selection - Round Two
ts_r2 = randperm(p/2);
winners = zeros(p/4, n);
for k = 2:2:p/2
    if tmp_fitness(ts_r2(k-1)) > tmp_fitness(ts_r2(k))
        winners(k/2, :) = winners_r1(ts_r2(k), :);
    else
        winners(k/2, :) = winners_r1(ts_r2(k-1), :);
    end
end
new_pop(1:p/4, :) = winners;
new_pop(p/2+1:3*p/4, :) = winners;

% Crossover
crossover = randperm(p/2);
children = zeros(p/4, n);
for k = 2:2:p/2
    parent1 = winners_r1(crossover(k-1), :);
    child = winners_r1(crossover(k), :);
    ndx = ceil(n*sort(rand(1, 2)));
    while ndx(1) == ndx(2)
        ndx = ceil(n*sort(rand(1, 2)));
    end
    tmp = parent1(ndx(1):ndx(2));
    for kk = 1:length(tmp)
        child(child == tmp(kk)) = 0;
    end
    child = [child(1:ndx(1)) tmp child(ndx(1)+1:n)];
    child(child == 0) = [];
    children(k/2, :) = child;
end
new_pop(p/4+1:p/2, :) = children;
new_pop(3*p/4+1:p, :) = children;

% Mutate
mutate = randperm(p/2);
num_mutate = round(mutate_rate*p/2);
for k = 1:num_mutate
    ndx = ceil(n*sort(rand(1, 2)));
    while ndx(1) == ndx(2)
        ndx = ceil(n*sort(rand(1, 2)));
    end
    if rand < 0.75 % swap segment between two custom
        new_pop(p/2+mutate(k), ndx(1):ndx(2)) = ...
            fliplr(new_pop(p/2+mutate(k), ndx(1):ndx(2)));
    else % swap two custom
        new_pop(p/2+mutate(k), [ndx(1) ndx(2)]) = ...
            new_pop(p/2+mutate(k), [ndx(2) ndx(1)]);
    end
end


