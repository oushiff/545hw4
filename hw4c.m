alpha_z = 25;
beta_z = 6;
alpha_x = 8;
y0 = 0;
x0 = 1;
z0 = 0;
N = 10;
g = 1;

c = [1.0000; 0.6294; 0.3962; 0.2494; 0.1569; 0.0988; 0.0622; 0.0391; 0.0246; 0.0155];
rho_sq = [ 41.6667; 16.3934; 6.5359; 2.5840; 1.0235; 0.4054; 0.1606; 0.0636; 0.0252; 0.0252]/1000;
w = [ 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];

dt = 0.001;





function dx = get_dx(x) 
    dx = - alpha_x  * x
end

function dy = get_dy(z)
    dy = z
end 

function dz = get_dz(y, z, f)
    dz = alpha_z * (beta_z * (g - y) - z) + f
end

function psi = get_psi(x, c)
    psi = exp( - ((x - c).^2) / (2 * rho_sq) )
end

function phi = get_phi(x, c)
    psi = get_psi(x, c)
    phi = psi * x * (g - y0) / sum(psi)
end

function f = get_f()
    f = phi.' * w
end
