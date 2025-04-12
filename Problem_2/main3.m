clear;
clc;
close all

R_earth = 6371302;
mu = 398600.4415 * 10.^9; % гравитационный параметр Земли
e = 0; % эксцентриситет
omega = 0; % аргумент перицентра

params = struct('j2_coef', 1082.8*10^(-6), 'earth_radius', R_earth);

%count = 172329; % количество итераций
%count = ceil(5619.36*46) + 1;
%count = (ceil(5619.36*46) + 1);
count = (ceil(120*60*60) + 1);
step = 1; % шаг в секундах


N_orb = 1; % количество плоскостей
N_sat = 3; % количество спутников на плоскости
inc = 98/180*pi; % наклонение
H = 460250.3338059597; % высота орбиты
phase = 0.8; % фаза

Omega0 = 0/180*pi; % долгота восходящего узла 1-ой орбиты
u0 = zeros(1, N_orb); % начальный аргумент широты 1-го спутника
u = ones(N_sat, N_orb);
N = N_sat * N_orb; % количество спуников

lat = zeros(count, N); % широты
lon = zeros(count, N);

for j = 2:N_orb
    u0(j) = u0(j-1) + phase * 2 * pi / N;
end

Omega = Omega0:pi/N_orb:(Omega0+pi*(N_orb-1)/N_orb);


for j = 1:N_orb
    u(1, j) = u0(j);
end

for j = 1:N_orb
    for k = 2:N_sat
        u(k, j) = u(k-1, j) + phase * 2*pi/N_sat;
    end
end


a = R_earth + H;
r0 = zeros(3, 1);
v0 = zeros(3, 1);

r = zeros(3, N, count);
v = zeros(3, N, count);
r_mean = zeros(3, N, count);
v_mean = zeros(3, N, count);

C = {'b','r','g','m','k'} ;
ci = 0;

for k_orb = 1:N_orb
    for k_sat = 1:N_sat
        r0(:, 1) = a * [cos(u(k_sat, k_orb))*cos(Omega(k_orb)) - sin(u(k_sat, k_orb))*cos(inc)*sin(Omega(k_orb));
                        cos(u(k_sat, k_orb))*sin(Omega(k_orb)) + sin(u(k_sat, k_orb))*cos(inc)*cos(Omega(k_orb));
                        sin(u(k_sat, k_orb))*sin(inc)];

        v0(:, 1) = sqrt(mu / a) * [- sin(u(k_sat, k_orb))*cos(Omega(k_orb)) - cos(u(k_sat, k_orb))*cos(inc)*sin(Omega(k_orb));
                                   - sin(u(k_sat, k_orb))*sin(Omega(k_orb)) + cos(u(k_sat, k_orb))*cos(inc)*cos(Omega(k_orb));
                                     cos(u(k_sat, k_orb))*sin(inc)];
        time = zeros(1, count);

        r(:, k_sat + N_sat * (k_orb - 1), 1) = r0;
        v(:, k_sat + N_sat * (k_orb - 1), 1) = v0;
        time(1) = 0;
        
        f = @func_J2;

        for j = 2:count

            [r(:, k_sat + N_sat * (k_orb - 1), j), v(:, k_sat + N_sat * (k_orb - 1), j)] =... 
                runge(f, mu, r(:, k_sat + N_sat * (k_orb - 1), j-1), v(:, k_sat + N_sat * (k_orb - 1), j-1), step, a);
            
            time(j) = time(j - 1) + step;
            [lat(j, k_sat + N_sat * (k_orb - 1)), lon(j, k_sat + N_sat * (k_orb - 1))] = ...
                 cart2geo(r(:, k_sat + N_sat * (k_orb - 1), j), time(j));
        end
     rr = squeeze(r(:, k_sat + N_sat * (k_orb - 1), :));

     ci = mod(ci + 1, 5) + 1;
     plot3(rr(1, :), rr(2, :), rr(3, :), 'Color', C{ci});
     hold on
     plot3(rr(1, 1), rr(2, 1), rr(3, 1), '.', 'MarkerSize', 20, 'Color', C{ci});
    end
end

hold off




% Визуализация треков на Земле
figure('Position', [100, 100, 1200, 600]);

% Вся Земля
subplot(1, 2, 1);
plot_earth();
hold on;
for idx = 1:N
    plot(lon(:, idx) * 180/pi, lat(:, idx) * 180/pi, '.', 'MarkerSize', 1);
end



sampling_step = 2809;
num_samples = floor(count / sampling_step);
sampled_lon = zeros(num_samples * N, 1);
index = 1;
% Заполняем массив выбранными значениями долготы
for sat = 1:N
    for j = 1:num_samples
        sample_index = j * sampling_step;
        if sample_index > count
            sample_index = count; % Чтобы не выйти за границы массива
        end
        sampled_lon(index) = mod(lon(sample_index, sat) + pi, 2*pi);
        index = index + 1;
    end
end


sorted_Lon = sort(sampled_lon );

distance_along_equator = sorted_Lon * (R_earth );

% 3. Находим разницы между соседними элементами
differences = diff(distance_along_equator);

%for idx = 1:31
%    plot(lon(558 * idx, 1) * 180/pi, lat(558 * idx, 1) * 180/pi, '.', 'MarkerSize', 10, "Color", "yellow");
%end

diff = abs(differences);
[max_diff, max_idx] = max(diff);
disp("max_diff");
disp(max_diff / 1000)

plot(lon(1, 1) * 180/pi, lat(1, 1) * 180/pi, '.', 'MarkerSize', 10, "Color", "green");
plot(lon(count, 1) * 180/pi, lat(count, 1) * 180/pi, '.', 'MarkerSize', 10, "Color", "magenta");
disp(lon(count, 1));

title('Треки спутников на всей Земле');
xlabel('Долгота (градусы)');
ylabel('Широта (градусы)');
grid on;

% Участок Земли [0;30] градусов долготы и [-20;20] градусов широты
subplot(1, 2, 2);
plot_earth();
hold on;
for idx = 1:N
    plot(lon(:, idx) * 180/pi, lat(:, idx) * 180/pi, '.', 'MarkerSize', 1);
end
xlim([0, 30]);
ylim([-20, 20]);
title('Треки спутников (долгота [0;30], широта [-20;20])');
xlabel('Долгота (градусы)');
ylabel('Широта (градусы)');
grid on;



save_interval = 60; 
save_steps = 1:save_interval:count;
num_saves = length(save_steps);

data = zeros(num_saves, 1 + 6*N);


for i = 1:num_saves
    idx = save_steps(i);
    data(i, 1) = (idx-1)*step; 
    
    for sat = 1:N
        col_offset = 1 + 6*(sat-1);
        data(i, col_offset+1:col_offset+3) = r(:, sat, idx)'; % Position
        data(i, col_offset+4:col_offset+6) = v(:, sat, idx)'; % Velocity
    end
end

% Create header
header = 't';
for sat = 1:N
    header = [header, sprintf(',x%d,y%d,z%d,vx%d,vy%d,vz%d', ...
              sat, sat, sat, sat, sat, sat)];
end

% Write to CSV file
filename = 'satellite_data_1.csv';
fid = fopen(filename, 'w');
fprintf(fid, '%s\n', header);
fclose(fid);

% Write data with 8 significant digits
dlmwrite(filename, data, '-append', 'precision', '%.8g', 'delimiter', ',');

disp(['Data saved to ' filename]);


initial_data = zeros(N, 6); % [x, y, z, vx, vy, vz] for each satellite
for sat = 1:N
    initial_data(sat, 1:3) = r(:, sat, 1)';
    initial_data(sat, 4:6) = v(:, sat, 1)';
    disp(r(:, sat, 1)');
    disp(v(:, sat, 1)');
end
disp(initial_data);
disp(93.565);
disp(46)


% Функция для визуализации Земли
function plot_earth()
    load coastlines; % Загрузка данных о береговой линии
    plot(coastlon, coastlat, 'k');
    hold on;
    axis equal;
    grid on;
end

function [lat, lon] = cart2geo(r, time)
    % Учитываем вращение Земли (угловая скорость в радианах/секунду)
    omega_earth = 7.292115e-5; 
    
    x = r(1);
    y = r(2);
    z = r(3);
    
    % Вычисляем "небесную" долготу (без учета вращения Земли)
    celestial_lon = atan2(y, x);
    
    % Корректируем долготу с учетом вращения Земли
    lon = celestial_lon - omega_earth * time;
    
    % Нормализуем долготу в диапазон [-pi, pi]
    lon = mod(lon + pi, 2*pi) - pi;
  %  lon = mod(lon, 2*pi);
    lat = atan2(z, sqrt(x^2 + y^2));
end
