LeakyIF_3();

dt = 0.0001;
input_currents.I_0=1e-8;

%psc = 0.08:0.09:dt;

%psc = 0:dt:dt;
%%psc = 0:dt:0.01;
psc = 0.05:dt:0.055;
%psc = 0.08:dt:0.6;
%psc = numpy.arange(0.02, dt+0.03, dt)

input_currents.psc=psc;
LeakyIF_3(input_currents,0.2,dt);