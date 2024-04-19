%% Passage Pipe Parameters (space between power and displacer pistons)
passage_pipe = struct();

% Geometry
passage_pipe.geom.radius = 4e-2; % m
passage_pipe.geom.length = 1e-2; % m 
passage_pipe.geom.area = pi*passage_pipe.geom.radius^2;
passage_pipe.geom.hydraulic_diam = 2*passage_pipe.geom.radius;


% Passage Pipe thermal_friction
passage_pipe.therm_fric.length_add = 0.01*passage_pipe.geom.length;
passage_pipe.therm_fric.roughness = 15e-6; %[m]
passage_pipe.therm_fric.Re_lam = 2000;
passage_pipe.therm_fric.Re_tur = 4000;
passage_pipe.therm_fric.shape_factor = 64;
passage_pipe.therm_fric.Nu_lam = 3.66;

%% Alternator Paramters
alternator = struct();

% Rotational Electromechanical Converter Parameters
alternator.converter.efficiency = 0.6;

% Electrical Load Parameters
alternator.load.resistance = 100; % ohms

%% Rhombic drive parameters:
% Geometry (geom)
rhombic_drive.geom.gear_rad = 8e-2; % m
rhombic_drive.geom.gear_mass = 15.; % kg
rhombic_drive.geom.crank_inertia = 0.5*rhombic_drive.geom.gear_mass*rhombic_drive.geom.gear_rad^2; % kg*m^2

% Drive paramters
rhombic_drive.drive.disp_link_length = 5.9953e-2; % m
rhombic_drive.drive.power_link_length = 8.0005e-2; % m 
rhombic_drive.drive.disp_crank_radius = 0.7716e-2; % m
rhombic_drive.drive.power_crank_radius = 1.3236e-2; % m

% Rotational damping
tau_damp = 10; %[s] time for inertia speed to drop to 25% of the initial speed
rhombic_drive.rot_damp = log(4)*rhombic_drive.geom.crank_inertia/tau_damp; % N*m*s/rad

%% Displacer Piston Parameters
displacer_piston = struct();
% Geometry (geom)
displacer_piston.geom.radius = 5e-2; % m
displacer_piston.geom.cyl_inner_radius = 4.95e-2; % m
displacer_piston.geom.length = 12e-2; % m
displacer_piston.geom.wall_thickness = 0.6e-2; % m
displacer_piston.geom.clearance_vol = 4e-6; % m^3

% Translational mechanical converter (piston) parameters
displacer_piston.piston.interf_area = pi*displacer_piston.geom.cyl_inner_radius^2; % m^2
displacer_piston.piston.ini_disp_cooler = rhombic_drive.drive.disp_crank_radius; % m
displacer_piston.piston.ini_disp_heater = -rhombic_drive.drive.disp_crank_radius; % m
displacer_piston.piston.dead_vol = displacer_piston.geom.clearance_vol + pi*(displacer_piston.geom.radius^2 - displacer_piston.geom.cyl_inner_radius^2)*displacer_piston.geom.length; % m^3
displacer_piston.piston.portA_area = -pi*(displacer_piston.geom.cyl_inner_radius^2 - displacer_piston.geom.radius^2); % m^2

% Translational hard stop (hardstop) parameters
displacer_piston.hardstop.upbound = 1.03*2*rhombic_drive.drive.disp_crank_radius; % m
displacer_piston.hardstop.lowbound = -1.03*2*rhombic_drive.drive.disp_crank_radius; % m
displacer_piston.hardstop.trans_region = 0.03*2*rhombic_drive.drive.disp_crank_radius; % m

%% Power Piston Parameters
power_piston = struct();
% Geometry (geom)
power_piston.geom.radius = 5e-2; % m
power_piston.geom.length = 12e-2; % m
power_piston.geom.clearance_vol = 4e-6; % m^3

% Translational mechanical converter (piston) parameters
power_piston.piston.interf_area = displacer_piston.piston.interf_area; % m^2
power_piston.piston.ini_disp = 0; % m
power_piston.piston.dead_vol = power_piston.geom.clearance_vol + pi*(power_piston.geom.radius^2 - displacer_piston.geom.cyl_inner_radius^2)*power_piston.geom.length; % m^3

% Translational hard stop (hardstop) parameters
power_piston.hardstop.upbound = 1.02*2*rhombic_drive.drive.power_crank_radius; %m
power_piston.hardstop.lowbound = -0.02*rhombic_drive.drive.power_crank_radius; % m
power_piston.hardstop.trans_region = 0.015*2*rhombic_drive.drive.power_crank_radius; % m


%% Ambient
geometry.cooler.length = 20e-2; % m

rad_ext = displacer_piston.geom.radius + displacer_piston.geom.wall_thickness;
rad_int = displacer_piston.geom.radius;
fins_equiv_area_factor = 4;

ambient.temperature = 293; % K
ambient.pext = 1; % MPa

ambient.fins.cp = 425; %J/kg/K Specific heat fins thermal mass
rho = 7000; % kg/m^3 steel density
finlength = 0.6*rad_ext;
ambient.fins.m = 0.35*rho*geometry.cooler.length*pi*((rad_ext+finlength)^2 - rad_ext^2); % kg mass fins

a.ambient.fins.h = 30; %W/(m^2*K)  Convection coefficient
a.ambient.fins.area = 2*pi*(rad_ext+finlength)*geometry.cooler.length*fins_equiv_area_factor; % m^2 Surface area
a.fins.gas.h = 60; % W/(m^2*K)  Convection coefficient
a.fins.gas.area = 2*pi*rad_int*geometry.cooler.length + pi*rad_int^2; % m^2 Surface area

clear rad_ext rad_int fins_equiv_area_factor rho finlength;

%% Initial state

state_init.T0 = 300; % K
state_init.p0 = 0.4298*ambient.pext; % MPa 


%% Flame

rad_ext = displacer_piston.geom.radius + displacer_piston.geom.wall_thickness;
rad_int = displacer_piston.geom.radius;
L = 40e-2;

flame.temperature = 2000; % K

flame.glass.cp = 700; % J/kg/K Specific heat glass thermal mass
rho = 4000; % kg/m^3 glass
flame.glass.m = 0.5*rho*L*pi*(rad_ext^2 - rad_int^2); % kg mass glass

f.flame.glass.area = 2*pi*rad_ext*(L - geometry.cooler.length) + pi*rad_ext^2; % m^2 Surface area
f.glass.gas.h = 60; % W/(m^2*K) Convection coefficient
f.glass.gas.area = 2*pi*rad_int*(L - geometry.cooler.length) + pi*rad_int^2; % m^2 Surface area

clear L rho rad_ext rad_int;

%Thermal masses temperature init
ambient.fins.T.init = 0.98*state_init.T0+0.02*flame.temperature; % Close to steady-state
flame.glass.T.init = 0.9*flame.temperature+0.1*ambient.temperature; % Close to steady-state


%% Regenerator Parameters
A = pi*(displacer_piston.geom.radius^2 - displacer_piston.geom.cyl_inner_radius^2); % Area m^2
P = 2*pi*(displacer_piston.geom.radius + displacer_piston.geom.cyl_inner_radius); % Perimeter m

% Geometry (geom)
regenerator.geom.length = displacer_piston.geom.length - 2*rhombic_drive.drive.disp_crank_radius ; %m
regenerator.geom.area = A; % m^2
regenerator.geom.hydraulic_diam = 4*A/P;

% Values (Re = Reynolds Number - Laminar Flow - Turbulent Flow) 
regenerator.therm_fric.length_add = 0.01*regenerator.geom.length;
regenerator.therm_fric.roughness = 15e-6; % m
regenerator.therm_fric.Re_lam = 2000;
regenerator.therm_fric.Re_tur = 4000;
regenerator.therm_fric.shape_factor = 64;
regenerator.therm_fric.Nu_lam = 3.66;
rho = 4000; % kg/m^3 glass

% External and Internal radius
rad_ext = displacer_piston.geom.radius + displacer_piston.geom.wall_thickness; % (Glass Thickness)
rad_int = displacer_piston.geom.radius;
% Mass
regenerator.m = 0.5*rho*pi*(rad_ext^2 - rad_int^2)*regenerator.geom.length;

% Conduction Cooler Part
regenerator.Cooler.area = regenerator.geom.area; 
regenerator.Cooler.thickness = regenerator.geom.length/2;
regenerator.Cooler.k = 80; % W/(m*K)

% Conduction Heater Part
regenerator.Heater.area = regenerator.geom.area;
regenerator.Heater.thickness = regenerator.geom.length/2;
regenerator.Heater.k = 80; % W/(m*K) %% Thermal Conductivity

regenerator.Tinit = 0.6*flame.glass.T.init+0.4*state_init.T0;  % Close to steady-state
clear P A rho rad_ext rad_int;

%% Impulse torque

dt = 1e-1;     % impulse duration (s)
deltaOm = 100;  % increase in angular speed (rad/s)

impulse_torque.t_start = 1; % s
impulse_torque.t_end = 1 + dt; % s
impulse_torque.torque = deltaOm*rhombic_drive.geom.crank_inertia/dt; % N*m

clear deltaOm dt;