%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HOUSEKEEPING & MATERIAL DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath("calfem/fem/"))  % CALFEM library

% --- Material properties (Table 1) ---
ECopper = 128e9;      % Copper E [Pa]
ENylon  = 3e9;        % Nylon E  [Pa]
vCopper = 0.36;       % Copper Poisson's ratio
vNylon  = 0.39;       % Nylon Poisson's ratio
alphaCopper = 17.6e-6; 
alphaNylon  = 80e-6;  % Thermal expansion [1/K]

% --- Convection and thermal properties ---
acopper   = 40;       % Convection coefficient for copper
rhoCopper = 8930;     % Density [kg/m^3]
rhoNylon  = 1100;     
cpCopper  = 386;      % Specific heat [J/kg-K]
cpNylon   = 1500;     
kCopper   = 385;      % Thermal conductivity [W/m-K]
kNylon    = 0.26;

% --- Subdomain property arrays ---
%  subdomain #4 => Nylon; the others => Copper
rho_subdomain    = [rhoCopper, rhoCopper, rhoCopper, rhoNylon,  rhoCopper, rhoCopper];
cp_subdomain     = [cpCopper,  cpCopper,  cpCopper,  cpNylon,   cpCopper,  cpCopper];
k_subdomain      = [kCopper,   kCopper,   kCopper,   kNylon,    kCopper,   kCopper];
alpha_subdomain  = [alphaCopper, alphaCopper, alphaCopper, alphaNylon, alphaCopper, alphaCopper];
v_subdomain      = [vCopper, vCopper, vCopper, vNylon, vCopper, vCopper];
E_subdomain      = [ECopper, ECopper, ECopper, ENylon, ECopper, ECopper];

% --- Other global parameters ---
T_inf     = 18;     % Ambient temperature [C]
h         = -1e5;   % Heat flux (negative => outward)
thickness = 100;    % For plane strain, thickness is somewhat arbitrary

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MESH & GEOMETRY SETUP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mesh = load('4nmesh.mat');   % Load PDETool-generated mesh
p = mesh.p;  % Node coordinates stored in columns [x; y]
e = mesh.e;  % Edges  (rows: nodeIDs, row5 => segment label)
t = mesh.t;  % Triangles (rows: nodeIDs, row4 => subdomain label)

coord = p';  
enod  = t(1:3,:)';  % Each row => 3 node indices for one element
nelm  = size(enod,1);
nnod  = size(coord,1);

% --- DOF numbering ---
dof   = (1:nnod)';         % For thermal: 1 DOF per node
dof_S = [(1:nnod)' , (nnod+1:2*nnod)'];  % For elasticity: 2 DOFs per node

% --- Create topology/connectivity matrices for CALFEM ---
for ie = 1:nelm 
    % edof_S => expands each element with global DOF indices for x & y
    edof_S(ie,:) = [ ie , ...
        dof_S(enod(ie,1),:), dof_S(enod(ie,2),:), dof_S(enod(ie,3),:) ];
    % edof => expands each element with node indices
    edof(ie,:)   = [ ie , enod(ie,:) ];
end

% --- Extract node coordinates per element [Ex, Ey] ---
[Ex, Ey] = coordxtr(edof, coord, dof, 3);

% --- Identify boundary edges from PDETool labeling ---
er = e([1 2 5],:);  % Keep rows 1,2 => nodeIDs, row5 => boundary label

% Segment labels from PDETool:
conv_segments           = [1 2 3 5 13 15 20 25]; 
flux_segments           = [18];
xdisplacement_segments  = [4 16 17 18];
ydisplacement_segments  = [16 17 18 21];

% Arrays to store node pairs belonging to each boundary condition
edges_conv          = [];    
edges_flux          = [];
edges_xdisplacement = [];
edges_ydisplacement = [];

% --- Translate PDETool segments => node pairs (edges) ---
for i = 1:size(er,2)
    seg_label = er(3,i);
    if ismember(seg_label, conv_segments)
        edges_conv = [edges_conv er(1:2,i)];
    end
    if ismember(seg_label, flux_segments)
        edges_flux = [edges_flux er(1:2,i)];
    end
    if ismember(seg_label, xdisplacement_segments)
        edges_xdisplacement = [edges_xdisplacement er(1:2,i)];
    end
    if ismember(seg_label, ydisplacement_segments)
        edges_ydisplacement = [edges_ydisplacement er(1:2,i)];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TASK (a) - STATIONARY HEAT CONVECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build conduction system (K,f) for conduction+convection+flux
[K, f] = buildConductionSystem(nnod, edof, Ex, Ey, coord, ...
    t, k_subdomain, edges_conv, edges_flux, acopper, T_inf, thickness, h);

% Solve stationary conduction
T = solveq(K, f);

% --- Plot stationary temperature ---
eT = extract(edof, T);
figure()
patch(Ex', Ey', eT')  % Quick patched surface
figure()
patch(Ex', Ey', eT','EdgeColor','none')
title('Temperature distribution [C]')
colormap(hot); 
cb = colorbar;
ylabel(cb,'°C','FontSize', 14);
xlabel('x-position [m]')
ylabel('y-position [m]')
axis equal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TASK (b) - TRANSIENT HEAT CONVECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) Build capacity matrix C for smaller time step (fine resolution)
C = buildCapacityMatrix(nnod, edof, Ex, Ey, thickness, ...
    rho_subdomain, cp_subdomain, t);

% 2) Time integration parameters
T_max    = 0.9 * max(T);  % 90% of max from stationary
time_step= 0.05;
end_time = 100;
nbr_steps= end_time/time_step + 1;

T_0 = 18;  % Initial temperature
T_transient = zeros(nnod, nbr_steps);
T_transient(:,1) = T_0;

% --- Forward Euler-like or trapezoidal step (depending on chosen approach) ---
for i = 2:nbr_steps
    T_previous = T_transient(:, i-1);
    T_transient(:,i) = (C + time_step*K) \ (time_step*f + C*T_previous);
end

% Identify time to reach 90% of max stationary
T_transient_max = max(T_transient);
time90  = 0;
index90 = 0;
for i = 1:nbr_steps
    if (T_transient_max(i) >= T_max)
        time90  = i * time_step;
        index90 = i;
        break
    end
end

% Plot snapshots at selected times
index3    = 0.03 * index90;  
snapshots = [1, ...
             round((1+index3/2)/2), ...
             round(index3/2), ...
             round((index3+index3/2)/2), ...
             floor(index3)];

for i = 1:5
    time = (snapshots(i)-1) * time_step;
    eT   = extract(edof, T_transient(:, snapshots(i)));
    figure()
    patch(Ex', Ey', eT','EdgeColor','none')
    title(['Temperature distribution [C], t= ', num2str(time), ' s'])
    colormap(hot);
    cb = colorbar;
    caxis([15, 28])
    ylabel(cb,'°C','FontSize',14);
    xlabel('x-position [m]')
    ylabel('y-position [m]')
    axis equal
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TASK (b) - TRANSIENT HEAT CONVECTION VIDEO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rebuild C for a coarser time step => for faster animation
C = buildCapacityMatrix(nnod, edof, Ex, Ey, thickness, ...
    rho_subdomain, cp_subdomain, t);

T_max     = 0.9 * max(T);
time_step = 0.5;
end_time  = 100;
nbr_steps = end_time/time_step + 1;
T_0       = 18;

T_transient = zeros(nnod, nbr_steps);
T_transient(:,1) = T_0;

figure()
eT = extract(edof, 18.*ones(nnod,1));
h  = patch('XData',Ex','YData',Ey','CData',eT','FaceColor','interp');
title('Time: 0','HorizontalAlignment','center','FontName','courier');
colormap(hot); 
colorbar; 
caxis([100 143.8]);

for i = 2:nbr_steps
    T_previous = T_transient(:,i-1);
    T_transient(:,i) = (C + time_step*K) \ (time_step*f + C*T_previous);
    eT = extract(edof, T_transient(:,i));

    set(h,'XData',Ex','YData',Ey','CData',eT');
    title(['Time (s): ', num2str((i-1)*time_step)], ...
           'HorizontalAlignment','center','FontName','courier');
    drawnow;

    if (max(T_transient(:,i)) > 141)
        pause(0.02)
    end
end

T_transient_max = max(T_transient);
time90 = 0;
index  = 0;
for i = 1:nbr_steps
    if (T_transient_max(i) >= T_max)
        time90 = i*time_step;
        index  = i;
        break
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TASK (c) - THERMOELASTICITY (TRANSIENT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ndof        = 2*nnod; 
u_transient = zeros(ndof, size(T_transient,2));

% --- Loop over time steps to compute displacements ---
for i = 1:size(T_transient,2)

    % Initialize global stiffness & force for elasticity
    K_el = zeros(ndof, ndof);
    f_el = zeros(ndof, 1);

    % Current temperature field
    T_current = T_transient(:, i);
    deltaT    = T_current - T_0; 

    % --- Build elasticity system (K_el, f_el) ---
    for elementnbr = 1:nelm
        subdomain     = t(4, elementnbr);  % Which subdomain this element is in
        current_E     = E_subdomain(subdomain);
        current_v     = v_subdomain(subdomain);
        current_alpha = alpha_subdomain(subdomain);

        D_constant = current_E / ((1+current_v)*(1-2*current_v));
        current_D  = D_constant .* [1 - current_v,  current_v,         0;
                                    current_v,      1 - current_v,     0;
                                    0,              0,        0.5*(1-2*current_v)];

        % Average element temperature => mean_deltaT
        mean_deltaT = (deltaT(t(1,elementnbr)) + ...
                       deltaT(t(2,elementnbr)) + ...
                       deltaT(t(3,elementnbr)) ) / 3.0;

        % Thermal stress vector
        D_epsilon_deltaT = current_alpha * mean_deltaT * current_D * [1;1;0];

        % Element matrices
        element_x = Ex(elementnbr,:);
        element_y = Ey(elementnbr,:);

        Ke       = plante(element_x, element_y, [2 thickness], current_D);
        f_deltaT = plantf(element_x, element_y, [2, thickness], D_epsilon_deltaT');

        % Assemble into global K_el, f_el
        indx = edof_S(elementnbr,2:end);
        K_el(indx,indx) = K_el(indx,indx) + Ke;
        f_el(indx)      = f_el(indx)      + f_deltaT;
    end

    % --- Apply displacement BC (edges_xdisplacement => u_x=0, edges_ydisplacement => u_y=0) ---
    bc = applyBoundaryDisplacements(edges_xdisplacement, edges_ydisplacement, nnod);

    % --- Solve for displacements at this time step ---
    [u, ~] = solveq(K_el, f_el, bc);
    u_transient(:, i) = u;
end

% --- Compute final-step von Mises stress ---
deltaT = T_transient(:, end) - T_0;
ed     = extract(edof_S, u_transient(:, end));  % Displacements in last time step

von_mises = zeros(1, nelm);
for elementnbr = 1:nelm
    subdomain     = t(4, elementnbr);
    current_E     = E_subdomain(subdomain);
    current_v     = v_subdomain(subdomain);
    current_alpha = alpha_subdomain(subdomain);

    D_constant = current_E / ((1+current_v)*(1-2*current_v));
    current_D  = D_constant .* [1 - current_v,  current_v,    0;
                                current_v,      1 - current_v,0;
                                0,              0,       0.5*(1-2*current_v)];

    mean_deltaT = (deltaT(t(1,elementnbr)) + ...
                   deltaT(t(2,elementnbr)) + ...
                   deltaT(t(3,elementnbr)) ) / 3.0;

    element_x  = Ex(elementnbr,:);
    element_y  = Ey(elementnbr,:);
    disp_local = ed(elementnbr,:);

    [es, ~] = plants(element_x, element_y, [2, thickness], current_D, disp_local);

    % Subtract pure thermal-stress part
    D_epsilon_deltaT = current_alpha * mean_deltaT * current_D * [1;1;0];
    es = es - D_epsilon_deltaT';  

    sigma_xx = es(1);
    sigma_yy = es(2);
    tau_xy   = es(3);

    % For plane strain => sigma_zz included
    sigma_zz = current_v*(sigma_xx + sigma_yy) - current_alpha*mean_deltaT*current_E;

    % Von Mises formula
    von_mises(elementnbr) = sqrt( sigma_xx^2 + sigma_yy^2 + sigma_zz^2 ...
                                - sigma_xx*sigma_yy - sigma_xx*sigma_zz ...
                                - sigma_yy*sigma_zz + 3*(tau_xy^2) );
end

% Compute nodal von Mises by averaging over connected elements
von_mises_nodes = zeros(nnod,1);
for n = 1:nnod
    [elemIDs,~] = find(edof(:,2:4) == n);
    von_mises_nodes(n) = sum(von_mises(elemIDs)) / numel(elemIDs);
end

% --- Plot final-step stress field ---
stress = extract(edof, von_mises_nodes);
figure()
patch(Ex', Ey', stress')
cb = colorbar;
ylabel(cb,'$\mathrm{N/m^2}$','Interpreter','latex','FontSize',14);
xlabel('$x\,(m)$','Interpreter','latex','FontSize',14)
ylabel('$y\,(m)$','Interpreter','latex','FontSize',14)
title('Effective von Mises stress field (final step)')
hold on

% Indicate max stress node
node_index = find(von_mises_nodes == max(von_mises_nodes));
pl = plot(coord(node_index,1), coord(node_index,2), 'x','MarkerSize',14, ...
          'LineWidth',2,'Color','b');
legend(pl,'Max von Mises stress','Location','best')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RENDERING: TEMPERATURE & DISPLACEMENT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Renderer', 'painters', 'Position', [100 100 900 600])
eT = extract(edof, 18.*ones(nnod,1));   % Just a baseline temperature for initial plot
patch(Ex',Ey',[0 0 0],'EdgeColor','none','FaceAlpha',0.3)
axis equal
title('Displacement field [Magnitude enhancement 15]')
h = patch('Xdata',Ex','YData',Ey','CData',eT','FaceColor','interp');

colormap(hot);
colorbar;
caxis([0 175]);

for i = 2:length(T_transient(1,:))   % Animate over time steps
    eT = extract(edof, T_transient(:,i));   % Temperature distribution at time i
    ed = extract(edof_S, u_transient(:,i)); % Displacements at time i
    
    % Magnify displacements for visualization
    mag  = 15;
    exd  = Ex + mag*ed(:,1:2:end);
    eyd  = Ey + mag*ed(:,2:2:end);

    set(h,'Xdata',exd','YData',eyd','CData', eT');
    drawnow;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOCAL HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [K, f] = buildConductionSystem(nnod, edof, Ex, Ey, coord, ...
    t, k_subdomain, edges_conv, edges_flux, acopper, T_inf, thickness, h)
% buildConductionSystem:
%   Builds the global stiffness matrix (K) and force (f) 
%   for the stationary conduction problem, accounting for:
%       - conduction within each subdomain
%       - convection on edges in 'edges_conv'
%       - flux boundary conditions on edges in 'edges_flux'
%
% Input:
%   - nnod, edof, Ex, Ey, coord: mesh data
%   - t, k_subdomain: subdomain info
%   - edges_conv, edges_flux: boundary edges from PDETool
%   - acopper, T_inf, thickness, h: physical/BC parameters
% Output:
%   - K, f: global conduction system

D = eye(2);
K = zeros(nnod);
f = zeros(nnod,1);

nelm = size(edof,1);

for elementnbr = 1:nelm
    % Identify subdomain => pick correct conductivity
    subdomain = t(4, elementnbr);
    current_k = k_subdomain(subdomain);
    current_D = current_k * D;  % [k 0; 0 k]

    % Basic conduction matrix via flw2te
    element_x = Ex(elementnbr,:);
    element_y = Ey(elementnbr,:);
    [Ke, fe]  = flw2te(element_x, element_y, thickness, current_D, 0);

    % Add convection if edges in 'edges_conv'
    for iConv = 1:length(edges_conv(1,:))
        if ismember([edges_conv(1,iConv), edges_conv(2,iConv)], edof(elementnbr,2:4))
            x1 = coord(edges_conv(1,iConv),1);
            y1 = coord(edges_conv(1,iConv),2);
            x2 = coord(edges_conv(2,iConv),1);
            y2 = coord(edges_conv(2,iConv),2);
            L  = sqrt((x1 - x2)^2 + (y1 - y2)^2);

            % Convection element matrix
            Kc = thickness*acopper*(L/6)*[0 0 0; 0 2 1; 0 1 2];
            fc = thickness*acopper*T_inf*(L/2)*[0 1 1]';
            Ke = Ke + Kc;
            fe = fe + fc;
        end
    end

    % Add flux if edges in 'edges_flux'
    for iFlux = 1:length(edges_flux(1,:))
        if ismember([edges_flux(1,iFlux), edges_flux(2,iFlux)], edof(elementnbr,2:4))
            x1 = coord(edges_flux(1,iFlux),1);
            y1 = coord(edges_flux(1,iFlux),2);
            x2 = coord(edges_flux(2,iFlux),1);
            y2 = coord(edges_flux(2,iFlux),2);
            L  = sqrt((x1 - x2)^2 + (y1 - y2)^2);

            % Negative sign (outward flux)
            fe = fe - thickness*h*(L/2)*[0 1 1]';
        end
    end

    % Assemble element into global
    indx = edof(elementnbr,2:end);
    K(indx,indx) = K(indx,indx) + Ke;
    f(indx)       = f(indx)       + fe;
end
end


function C = buildCapacityMatrix(nnod, edof, Ex, Ey, thickness, ...
    rho_subdomain, cp_subdomain, t)
% buildCapacityMatrix:
%   Builds the global capacity matrix C for transient conduction
%   using "plantml" from CALFEM.
%
% Input:
%   - nnod, edof, Ex, Ey: mesh data
%   - thickness: thickness (for plane problems)
%   - rho_subdomain, cp_subdomain: material density & heat capacity
%   - t: subdomain labels
% Output:
%   - C: global capacity matrix

C    = zeros(nnod, nnod);
nelm = size(edof,1);

for elementnbr = 1:nelm
    subdomain  = t(4, elementnbr);
    current_rho= rho_subdomain(subdomain);
    current_cp = cp_subdomain(subdomain);

    element_x = Ex(elementnbr,:);
    element_y = Ey(elementnbr,:);

    % plantml => capacity (mass * heat capacity) for each element
    Ce = plantml(element_x, element_y, current_rho*current_cp*thickness);

    % Assemble into global C
    indx = edof(elementnbr,2:end);
    C(indx,indx) = C(indx,indx) + Ce;
end
end


function bc = applyBoundaryDisplacements(edges_xdisp, edges_ydisp, nnod)
% applyBoundaryDisplacements:
%   Creates a boundary-condition (bc) array for zero displacement
%   along edges_xdisp => x-DOF = 0
%   along edges_ydisp => y-DOF = 0
%
% Output bc is [dof_id, value_of_bc]

bc = [ ...
    edges_xdisp(1,:)'; 
    edges_xdisp(2,:)';
    nnod + edges_ydisp(1,:)'; 
    nnod + edges_ydisp(2,:)' ...
    ];

bc = unique(bc);  % remove duplicates
bc = [bc, zeros(size(bc))];
end
