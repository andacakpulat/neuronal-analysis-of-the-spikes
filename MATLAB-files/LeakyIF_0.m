function [time, U_plot]=LeakyIF_0(input_currents,duration,dt)
%_____________________________________________________________________
%Note before start: This is a function. Translated to Python will be "def
%LeakyIF_0(input_currents,duration,dt):" and  indenting the code below.
%"[time, U_plot]" are two returned arrays each (see meaning below), Matlab stores all function
% outputs in this 1-D array (like in Numpy, or a list).
% Matlab does not need for us to "return" them in the code below. In Python, we would need 
% to return them. 
%_____________________________________________________________________
%
%Leaky integrated and fire model for a single neuron-only membrane
%potential equation
%Based on Prof. Kevin Gurney material.
%Adapted by Emili B-B, Bournemouth University. 
%Neuronal Analytics Masters Course.
%Suboptimal code, desgined for education/clarity/portability to Python or other languages,
%not for performance.
%
%close all;clc;
%Optional inputs: 
%           input_currents.psc=1-D array (vector), in seconds. Each entry
%           indicates when an excitatory current enters the cell (values from 0 to
%          the spike train duration, in intervals of dt). 'psc' stands for
%          'post-synaptic current', that is the current delivered from
%          another cell that spiked in the past. Note that
%          "input_currents" is an object which has two fields, this one
%          ("psc") and the one below ("I_0").
%
%          input_currents.I_0= Indicates the input currrent value. Note that the membrane potential increase caused by this current,
%          is just the intensity times the membrane resitance R, which is fixed here (Sessions I.3.1-I.3.2 slides). Suggest around 1e-8 Amp (that
%           would render an increase in membrane potential of 0.3 V per
%           incoming excitatory spike, or, 300 mV. given the membrane
%           resistance value hard-coded here, see below these comments).
%
%
%           duration=spike train duration i.e., duration of the experiment
%           (in seconds). The times in which a spike will produce an excitatory input
%           to the neuron are stored in a field "psc" of the object
%           "input_currents", that is,
%           "input_currents.psc".
%Outputs:
%            time=Vector of times in sec., e.g., for plotting outside this
%            function.
%            U_plot=Membrane potential values (that is, the output of the neuron showing spikes if any), of
%            size 1 x length(time).
%Note that in matlab there is no need to specify the variables to return
%in the function body (just in the definition).
%
% Last modification: 11/11/2022.

%% PARAMETER SETUP

%1. Fixed parameters: Membrane constants (all IS units). These indicate the cell
%morphology and the maximum temporal resolution of the sinmulation, and thus are hardcoded.
    
    %1.1 Basic parameters
    tau = 0.020; %Membrane time constant in seconds (R x C in the model shown in the slides)
    R = 3e7;%Membrane resistence in ohms
    U_rest = -0.07;% Resting potential of the cell in volts
    %
    %1.2 Simulation parameters
    if nargin<3
        dt = 0.0001;%0.0001; %Time step in seconds, i.e. temporal resolution of the integration.
        %Should be naturally very small, much less than tau. Suggest not to
        %modify it.
    end
    
    if nargin<2
        T = 0.2;%Duration of the period we want to simulate in seconds
    else
        T=duration;
    end
    
%2. Parameters to modify: 
    %These are just hard-coded examples of input pulses charging the membrane
    %capacitor to demostrate spiking
    % Pulses can be provided externally (see Optional inputs in this file) or hardcoded here
     if nargin<1    
         U_0=0.3;%In volts
         I_0=U_0/R;%1e-8 Amp.  
         %An example of an inital input for a vert short period
         psc =0:dt:dt; %This is like a Numpy 1-D array: NumPy.arange(0., dt, dt)
     else
         psc=input_currents.psc;%In Python, to copy the array but not referencing 
         %                                    it would be
         %                                    psc=NumPy.array(input_currents.psc, copy=true)
         %                                    or equivalent
         I_0=input_currents.I_0; %same, or you can use deepcopy() to copy the full object 
         %                                  "input_currents"
    end

%% MODEL

%1-Initializations (suggest do not modify)
    n_pcs = length(psc); %number of pre-synaptic excitatory currents (that is, positive input pulses delivered)
    index_pscs = round(psc ./ dt);%Indexes of the these events (0, 1, 2,... see below)
    n_steps = round(T ./ dt); %Number of steps in the simulation
    U = zeros(1, n_steps + 1); %Output voltage matrix
    U_plot = zeros(1, n_steps + 1); %Output voltage matrix -duplicated just for convenience in 
    %                                              plotting spikes, see below
    U(1) = U_rest; %The resting potential is the first value of membrane voltage
    U_plot(1)=U_rest;
    I = zeros(1, n_steps + 1);%Net input intenstiy value a each step of the simulation (see below)
    time = linspace(0, T, n_steps + 1);%Creates a vector of times for plotting. Same in python, Numpy.linspace...

%% SIMULATION
    %This is a "for" loop 
    for i=1:n_steps %We start from the beggining of the stimulation period T (i=1) and loop until the end of it 
        for k=1:n_pcs %Then, for each time step "i", loop over the times in which there is a pulse delivered, and check if it matches the current time
            if i == index_pscs(k)%We just loop this vector containing the time steps in which a current should be delivered
                %when the index in this vector corresponds to the actual
                %time step "i" (that is, the current time), a burst of
                %current of value I_0 is provided. This implementation is
                %not efficent and is done in this way just for clarity.
                %Note that in matlab the first index is "1", unlike
                %in Python. This is probably cleaner with a dictonary etc. if
                %you do not want to use nested loops, but it does the job
                %as it is now
                I(i) = I(i) + I_0;%Excitatory current burst -remember, in python is I[i] +=I_0, and the same with all parenthesis below
            end
        end
        dU =(dt/tau).*(U_rest-U(i)+I(i).*R);%Increment of membrane potential in the passive membrane equation
        %See I.3.1 last slides in
        %"Model-Mechanisms_NA_I_Intro_Neuronal_....pptx" file (remember in
        %Python in U[i] etc. instead)
        U(i+1) = U(i) + dU; 
        U_plot(i+1)=U(i+1);%Having two identical separare variables seems pointless here, is just for plotting.
        %It will become clearer why in the more advanced model versions ('LeakyIF_1.m' etc.)
    end
    
    %% PLOTTING
    figure(1)
    plot(time, U_plot, 'color', 'black');
    title('Membrane equation model');
    xlabel('time (s)')
    ylabel('voltage (U)')
    grid on