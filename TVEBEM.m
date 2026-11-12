N0 = 12; M0 = 24; % number of elements in theta (N0) and phi (M0) directions
Q = N0*M0; 
ra = 0.94; % radius of sphere at which field points are taken

om = 10; % angular frequency

rad = 1; % radius of each sphere
d = 2.1*rad; % centre-to-centre separation between spheres
V1 = [0, -d/2, 0]; V2 = [0, d/2, 0]; % position of centres

elasticR = [20, 10, 15, 2]; % mechanical properties in the order specified below
K = elasticR(1); % bulk modulus
muR = elasticR(2); % long-term shear modulus
mu0 = elasticR(3); % short-term shear modulus
rho = elasticR(4); % mass density

thermal = [0.003, 2, 10, 300]; %thermal properties in the order specified below
al = thermal(1); % coefficient of thermal expansion
kth = thermal(2); % thermal conductivity
c = thermal(3); % specific heat
T0 = thermal(4); % initial temp

tr = 1; % relaxation time

%frequency dependent shear modulus
mu = muR - (mu0-muR)*(1i*om*tr)/(1-1i*om*tr);

% obtain wavenumbers
l = K - 2/3*mu; 
A = l + 2*mu;
C = al*K*T0/A;
k1 = sqrt(1i*om*c*rho/kth); k2 = sqrt(1i*om*al*K/kth);
k3 = sqrt(rho*om^2/A);
    
a = 1/2*(k1^2+k3^2+C*k2^2); b = sqrt(a^2-k1^2*k3^2);
kc = sqrt(a-b); % thermal dominated wavenumber
kt = sqrt(a+b); % compression dominated wavenumber
ks = sqrt(rho*om^2/mu); % shear wavenumber

B = al*K; 
elastic = [mu, l, rho];
thermal = [al, kth, c, T0];
wavenumbers = [k1, k2, k3, kt, kc, ks];

elmt = [0;0;0];
    
for void = -1:2:1
    for na = 0:N0-1
        for ma = 1:M0
            elmt = cat(2, elmt, [na; ma; void]);
        end
    end
end

mat = cell(length(elmt)-1);
vec = zeros(length(elmt)-1, 1);
tvec = zeros(length(elmt)-1, 1);
Q_ind = numel(1:length(elmt)-1);


parfor a = 1:Q_ind
    nma = elmt(:, a+1);
    
    [xa1, xa2, xa3] = y(0, 0, quad(nma(1), nma(2), ra, N0, M0, nma(3)*d/2));

    % incident field
    vec(a) = -1i*exp(1i*kt*xa3)/kt; % displacement 
    tvec(a) = 1/(C*kt^2)*(k3^2-kt^2)*exp(1i*kt*xa3); % change in temp

    for q = 1:Q_ind
        nmq = elmt(:, q+1);     
        mat{a, q} = H_mat(nma(1), nma(2), nmq(1), nmq(2), N0, M0, ra, elastic,thermal, wavenumbers, B, nma(3), nmq(3), d); 

    end
end

fin_vec = zeros(4*length(elmt)-4, 1);
fin_vec(3:4:end) = vec(:);
fin_vec(4:4:end) = tvec(:);
fin_mat = cell2mat(mat);

uq = fin_mat\ fin_vec;

writematrix(uq, "output.txt")



function mat = H_mat(na, ma, nq, mq,nmax, mmax, ra,elastic, thermal,wavenumbers,B, sign1, sign2, d)
    mat = zeros(3,4);

    for i = 1:3
        for j = 1:3
            mat(i,j) = h(i, j, na, ma, nq, mq,nmax, mmax, ra,elastic, wavenumbers,B, sign1, sign2, d);
        end
         
        mat(i, 4) = h2(i,na, ma, nq, mq,nmax, mmax, ra,elastic, wavenumbers,B, sign1, sign2, d);
    end
    
    for j = 1:4
        mat(4, j) = h3(j, na, ma, nq, mq,nmax, mmax, ra,elastic,thermal,B, wavenumbers, sign1, sign2, d);
    end
    
end


function f1 = h(i, j, na, ma, nq, mq,nmax, mmax, ra,elastic, wavenumbers,B, sign1, sign2, d)
    
    dif = 0.2;

    Eta = -1:dif:1; Xi = -1:dif:1;
    [eta, xi] = meshgrid(Eta, Xi);
    
    [Y1, Y2, Y3] = y(0, 0, quad(na, ma, ra, nmax, mmax,sign1*d/2)); % field points
    [X1, X2, X3] = y(xi, eta, quad(nq, mq, 1, nmax, mmax,sign2*d/2)); % source points

    integrand = Uij_vector_BIR(i, j, Y1, Y2, Y3, X1,X2,X3, elastic, wavenumbers, B, nmax,mmax, nq,mq, eta);
    f1 = trapz(Eta, trapz(Xi, integrand,2));
end

function f2 = h2(j,na, ma, nq, mq,nmax, mmax, ra,elastic,wavenumbers,B, sign1, sign2, d)

    
    dif = 0.2;
    
    Eta = -1:dif:1; Xi = -1:dif:1;
    [eta, xi] = meshgrid(Eta, Xi);
    
    [Y1, Y2, Y3] = y(0, 0, quad(na, ma, ra, nmax, mmax,sign1*d/2)); % field points
    [X1, X2, X3] = y(xi, eta, quad(nq, mq, 1, nmax, mmax,sign2*d/2)); % source points

    int2 = temp_vector_BIR(j, Y1, Y2, Y3, X1,X2,X3, elastic, B,wavenumbers, nmax,mmax, nq,mq, eta);
    f2 = trapz(Eta, trapz(Xi, int2,2));
end

function f1 = h3(i, na, ma, nq, mq,nmax, mmax, ra,elastic,thermal,B, wavenumbers, sign1, sign2, d)
    
    dif = 0.2;

    Eta = -1:dif:1; Xi = -1:dif:1;
    [eta, xi] = meshgrid(Eta, Xi);
    
    [Y1, Y2, Y3] = y(0, 0, quad(na, ma, ra, nmax, mmax,sign1*d/2)); % field points
    [X1, X2, X3] = y(xi, eta, quad(nq, mq, 1, nmax, mmax,sign2*d/2)); % source points

    integrand = scalar_rigid_BIR(i, Y1, Y2, Y3, X1,X2,X3, elastic,thermal,B, wavenumbers, nmax,mmax, nq,mq, eta);
    f1 = trapz(Eta, trapz(Xi, integrand,2));
end

function [x1, x2, x3] = y(xi, eta, y_quad) 
    y1 = y_quad(:, 1); y2 = y_quad(:, 2); y3 = y_quad(:,3);
    x1 = ytest(xi, eta, y1); 
    x2 = ytest(xi, eta, y2); 
    x3 = ytest(xi, eta, y3);
end

function x = ytest(xi, eta, yi)    
    x = 1/4.*(1-xi).*(1-eta).*yi(1) +  1/4.*(1+xi).*(1-eta).*yi(2) +  1/4.*(1-xi).*(1+eta).*yi(3) + 1/4.*(1+xi).*(1+eta).*yi(4);
end

function q = quad(n, m, r, nmax, mmax,shift)
    % q is the vector of nodes for element (n,m) on sphere x^2+(y-shift)^2 + z^2 = r^2
    % where n = 0, ..., nmax - 1 and m = 1, ..., mmax
    q = [node(n, m, r, nmax, mmax,shift); node(n, m+1, r, nmax, mmax,shift); node(n+1, m, r, nmax, mmax,shift); node(n+1, m+1, r, nmax, mmax,shift)];
end

function x = node(n, m, r, nmax, mmax,shift)
    % position of the (n,m)th node on sphere x^2+(y-shift)^2 + z^2 = r^2
    % where n = 0, ..., nmax - 1 and m = 1, ..., mmax
    x = [r.*sin(pi*n/nmax).*cos(2*pi*m/mmax), r.*sin(pi*n/nmax).*sin(2*pi*m/mmax)-shift, r.*cos(pi*n/nmax)];
end