% Generate x, y ordered pairs for y = sin(x)
x = linspace(0, 2*pi, 100);
y = sin(x);

% Salt the ordered pairs
salt_range = 0.2;
salted_y = y + salt_range * rand(size(y)) - salt_range/2;

% Smooth the ordered pairs using a moving average filter
window_size = 5;
smoothed_y = filter(ones(1, window_size)/window_size, 1, salted_y);

% Plot the original points
subplot(3, 1, 1);
plot(x, y, '-');
title('Original Points');
xlabel('x');
ylabel('y');

% Plot the salted points
subplot(3, 1, 2);
plot(x, salted_y, '-');
title(sprintf('Salted Points (range = %f)', salt_range));
xlabel('x');
ylabel('y');

% Plot the smoothed points
subplot(3, 1, 3);
plot(x, smoothed_y, '-');
title(sprintf('Smoothed Points (window size = %d)', window_size));
xlabel('x');
ylabel('y');
