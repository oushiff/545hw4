input = importdata('imitation.data');
Y = input(:, 1);
YD = input(:, 2);
YDD = input(:, 3);
len = length(X);

global alpha_z beta_z alpha_x g c rho_sq x0 y0;
y0 = 0;
x0 = 1;
z0 = 0;
alpha_z = 25;
beta_z = 6;
alpha_x = 8;
c = [1.0000; 0.6294; 0.3962; 0.2494; 0.1569; 0.0988; 0.0622; 0.0391; 0.0246; 0.0155];
rho_sq = [ 41.6667; 16.3934; 6.5359; 2.5840; 1.0235; 0.4054; 0.1606; 0.0636; 0.0252; 0.0252]/1000;
y0 = Y(1);
N = 10;
g = Y(1001);
t = 0;
td = 0.001;
W = [];
i = 2;
while i <= 1001
    x = get_x(t);
    phi = get_phi(x);
    f = YDD(i) - alpha_z * (beta_z * (g - Y(i)) - YD(i));
%     w =  inv(phi * (phi.')) * phi * f;
    w = f / phi ;
    W = [W ; w];
    i = i + 1;
    t = t + td;
end
w = (sum(W)/len).'
x = x0;
y = y0;
yd = z0;
ydd = 0;
z = z0;
t = 0;
new_Y = [];
new_YD = [];
new_PSI = [];
new_YDD = [];
while t <= 1
    new_Y = [new_Y; y];
    new_YD = [new_YD; yd];
    psi = get_psi(x);
    new_PSI = [new_PSI; psi.'];
    ydd = get_zd(y, z, get_f(x)); 
    new_YDD = [new_YDD; ydd];
    yd = yd + ydd * td;
    y = y + yd * td;
    z = yd;
    x = x + get_xd(x) * td;
    t = t + td;
end
figure
plot(0.001:0.001:1,[new_Y new_YD new_YDD])
%plot([new_Y new_YD new_YDD])
legend('y', 'yd', 'ydd')
figure
PSI;
plot(0.001:0.001:1,[new_PSI])
legend('psi1', 'psi2', 'psi3','psi4', 'psi5', 'psi6', 'psi7','psi8', 'psi9','psi10')


function x = get_x(t)
    global alpha_x x0;
    x = x0 * exp( - alpha_x * t);
end

function xd = get_xd(x) 
    global alpha_x;
    xd = - alpha_x  .* x;
end

function yd = get_yd(z)
    yd = z;
end 

function zd = get_zd(y, z, f)
    global alpha_z beta_z g;
    zd = alpha_z * (beta_z * (g - y) - z) + f;
end

function psi = get_psi(x)
    global c rho_sq;
    psi = exp( - ((x - c).^2) ./ (2 * rho_sq) );
end

function phi = get_phi(x)
    global g y0;
    psi = get_psi(x);
    phi = psi * x * (g - y0) / sum(psi);
end

function f = get_f(x)
    global w;
    f = get_phi(x).' * w;
end