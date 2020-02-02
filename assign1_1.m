% Jakob Horvath, u1092049
% Plots the Hubbard equation with different values of time, as well as the exponential function, e.

% show the hubbard equation evolving thru time, using four different 't'
% values
for t=0:0.25:0.75
    x = linspace(0, 1, 50);
    y = ones(50, 1) .* hubbard(x, t);
    plot(x, y), title('Hubbard w/ n = 50, t = 0.0, 0.25, 0.50, 0.75');    
    hold on;
end

% show a single hubbard equation's graphical output when 't' = 0.75
x = linspace(0, 1, 50);
y = ones(50, 1) .* hubbard(x, .75);
figure, plot(x, y), title('Hubbard w/ n = 50, t = 0.75');

% plot the exponential function, 'e'
x = linspace(0, 1, 100);
y = exp(x);
figure, plot(x, y), title('e');