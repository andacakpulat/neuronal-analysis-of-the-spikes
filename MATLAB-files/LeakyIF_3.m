function [time, U_plot]=LeakyIF_3(input_currents,duration,dt)
%_____________________________________________________________________
%Note before start: This is a function. Translated to Python will be "def
%LeakyIF_3(input_currents,duration,dt):" and  indenting the code below.
%"[time, U_plot]" are two returned arrays each (see meaning below), Matlab stores all function
% outputs in this 1-D array (like in Numpy, or a list).
% Matlab does not need for us to "return" them in the code below. In Python, we would need 
% to return them. 
%_____________________________________________________________________
%
%Leaky integrated and fire model for a single neuron-with
%refractory period, background noise and a primitive model of synapse:
%    1) excitatoy/inhibitory input currents (modelling AMPA-NMDA/GABA gates) and 
%    2) and adptive threshold (modelling depletion and recovery of synaptic resources) 
% 
%Based on Prof. Kevin Gurney material.
%Adapted by Emili B-B, Bournemouth University. 
%Neuronal Analytics Masters Course.
%Suboptimal code, desgined for education/clarity/portability to Python or other languages,
%not for performance.
%

%NOTE PLEASE: For differences with LeakyIF_2, search for the "Addition" word.

%close all;clc;
%Optional inputs: 
%           input_currents.psc=1-D array (vector), in seconds. Each entry
%           indicates when an excitatory current enters the cell (values from 0 to
%          the spike train duration, in intervals of dt). 'psc' stands for
%          'post-synaptic current', that is, the current delivered from
%          another cell that spiked in the past. Note that
%          "input_currents" is an object which has two fields, this one
%          ("psc") and the one below ("I_0").
%
%           Addition with respect to LeakyIF_2: 
%           input_currents.U_0=To reflect the idea that there is no true input current, membrane 
%                                       potential perturbation of 0.3 V
%                                       (300 mV).
%
%           duration=spike train duration i.e., duration of the experiment
%           (in seconds). Note that the times in which a spike will produce an excitatory input
%           to the neuron are stored in "input_currents.psc".
%
%           Addition with respect to LeakyIF_2:
%           input_currents.ipsc=Vector array, in seconds, each entry
%           indicates when an inhibitory current enters the cell (values from 0 to
%          the spike train duration, in intervals of dt). See I.3.2 slides.
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
%morphology and the maximum temporal resolution of the sinmulation, and thus are hardcoded

    %1.1 Basic 
    tau = 0.020; %Membrane time constant in seconds (R x C in the model)
    R = 3e7;%Membrane resistence in ohms
    U_rest = -0.07;% Resting potential of the cell in volts (also reset ponential)
    theta = -0.030;%Threshold of depolarizarion for a spike to trigger in volts
    %                    Typically 20-30 mv above  resting potential
    spikeVolt=0.1;%Spike value in volts. This is hard-coded since the I&F is too simple to reproduce a spike
    %
    %1.3 Simulation parameters
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
    
    %1.3 Refined parameters
      arp = 0.008;  %Absolute refactory period in seconds. Usally very small, << 10 ms
      backgroundI=3e-9;%Noise input backgound in Amp, present all the time from the background input in the medium
      
      %Additions with respect to LeakyIF_2: 
      %
      tau_adapt=0.15; %That is how fast the adaptive threshold decays to the inital " theta" value
      %when no spikes
      increase_threshold=0.012; %In volts, this is how much the threshold increses after a spike
      %This will make harder to initiate a spike after a previous spike+ the absolute refractory period window
      
%2. Parameters to modify: 
    %These are just examples of input pulses chargign the membrane capacitor to demostrate spiking
    % Pulses can be provided externally (see Optional inputs in this file) or hardcoded here
     ipsc =[]; %Addition with respect to LeakyIF_2: potential 'inhibitory' currents: input currents that generally do 
                  % not favour spiking even if large, tend to hyperpolarize membrane.
    if nargin<1    
         %Addition with respect to LeakyIF_2: Note that now, to be
         %consistent with the idea that the synapse is chemical and not
         %electrical, the input rather causes an alteration in the membrane potential, U_0. 
         U_0=0.3;%In volts
         I_0=U_0/R;%1e-8 Amp. 
         psc =0:dt:dt; %This emulates a failed initiation. 
         %See "LeakyIF_0" for more details and Python analogy
    else
         psc=input_currents.psc; %See "LeakyIF_0" for more details and Python analogy

         %Addition with respect to LeakyIF_2: Note that now, to be
         %consistent with the idea that the synapse is chemical and not
         %electrical, the input rather causes an alteration in the membrane potential, U_0
          I_0=input_currents.U_0/R; 

         %Addition with respect to LeakyIF_2: 'inhibitory' currents: input currents that generally do not favour spiking even if large, 
         %tend to hyperpolarize membrane
         ipsc=input_currents.ipsc; 
    end

%% MODEL

%1-Initializations (suggest do not modify)
    n_pcs = length(psc);%number of pre-synaptic excitatory currents (that is, positive input pulses delivered)
    n_ipcs = length(ipsc);%Addition with respect to LeakyIF_2: number of pre-synaptic inhibitory currents (that is, negative input pulses delivered)
    % index for each event
    index_pscs = round(psc ./ dt);%Indexes of the these events (0, 1, 2,... see below)
    index_ipscs = round(ipsc ./ dt);%Addition with respect to LeakyIF_2: same for inhibitory currents: prevent spiking.
    n_steps = round(T ./ dt); %Number of steps in the simulation
    U = zeros(1, n_steps + 1); %Output voltage matrix
    U_plot = zeros(1, n_steps + 1); %Output voltage matrix -duplicated just for convenience in 
    %                                              plotting spikes, see below.
    U(1) = U_rest; %The resting potential is the first value of membrane voltage
    U_plot(1)=U_rest;
    I = zeros(1, n_steps + 1);%Net input intenstiy value a each step of the simulation (see below)
    t_spike = 0; %Intial time since last spike
    n_spikes = 0; % Total number of spikes
    time = linspace(0, T, n_steps + 1);%Creates a vector of times for plotting
   randI = backgroundI* random('Normal', 0, 1, [1, n_steps]);%Background noise drawn form a normal distribution (central limit)
   %
   %Addition with respect to LeakyIF_2: Adapive threshold for spike increases after a spike and
    %decays to the default vale otherwise.
    theta_adapt= theta*ones(1, n_steps + 1);

%% SIMULATION

    for i=1:n_steps %We start from the beggining of the stimulation period T (i=1) and loop until the end of it 
        for k=1:n_pcs %Then, for each time step "i", loop over the times in which there is a pulse delivered, and check if it matches the current time
            if i == index_pscs(k)%See further details in LeakyIF_0.m
                I(i) = I(i) + I_0;%Excitatory current -remember, in python is I[i] +=I_0, and the same with all parenthesis below
            end
        end
        
        %Addition with respect to LeakyIF_2: Inhibitory current. See I.3.2
        %slides.
        for x=1:n_ipcs %Same for the inhibitory currents
            if i == index_ipscs(x)
                I(i) = I(i) - I_0;
            end
        end
        dU =(dt/tau).*(U_rest-U(i)+I(i)*R + randI(i)*R);%Passive membrane equation with small background input current
        U(i+1) = U(i) + dU;
        U_plot(i+1)=U(i+1);
        %Spike detection

        %Addition with respect to LeakyIF_2: adaptive threshold 
        if (U(i+1) > theta_adapt(i)) 
            if (n_spikes>0)
                if (time(i)>=(t_spike+arp)) 
                    U_plot(i+1)=spikeVolt;
                    U(i+1) = U_rest;%Reset membrane voltage to resting value.We could say here instead i+2 to be precise with the equation, but then we need to stop the i-th loop in n_steps-1
                    t_spike = time(i);
                    n_spikes = n_spikes+1;

                     %Addition with respect to LeakyIF_2: threhold increases after spike
                    theta_adapt(i+1)=theta_adapt(i)+increase_threshold;
                else
                    U(i+1) = U(i);
                    U_plot(i+1)=U(i+1);

                    %Addition with respect to LeakyIF_2: The threshold for spike decays to the
                    %default vale if no spike. Note that is the same equation as the passive membrane equation for U when no input.
                     theta_adapt(i+1) = theta_adapt(i)+(dt/tau_adapt)*(theta-theta_adapt(i));            
                end   
                    
            else
                U_plot(i+1)=spikeVolt;
                U(i+1) = U_rest;
                t_spike = time(i);
                n_spikes = n_spikes+1;
                %Addition with respect to LeakyIF_2: threshold increases after spike
                theta_adapt(i+1)=theta_adapt(i)+increase_threshold;
            end    
        else
                %Addition with respect to LeakyIF_2: The threshold for spike decays to the
                %default vale if no spike. Is the same equation for U
                %when no input.
                 theta_adapt(i+1) = theta_adapt(i)+(dt/tau_adapt).*(theta-theta_adapt(i));           
        end
         
        
    end
    
    %% PLOTTING
    figure
    plot(time, U_plot, 'color', 'black');
    hold on
    line([0,time(end)],[theta,theta],'Color','red','LineStyle','--');
    title_text_1='Leaky I&F with background current, refractory period, excit/inhib synapses';
    title_text_2=[' and adaptive threshold. Rest value: ',num2str(theta*1000),' mV. Refractory period ',num2str(1000*arp),' ms'];
    title({title_text_1,title_text_2});
    
    xlabel('time (s)')
    ylabel(['voltage (U), number of spikes=',num2str(n_spikes)])
    grid on
    
    %%Addition with respect to LeakyIF_2: plot the threshold for spiking
    %%dynamics
    
    figure
    plot(time, theta_adapt, 'color', 'red');
    hold on
   line([0,time(end)],[theta,theta],'Color','red','LineStyle','--');
    title(['Rest value of threshold: ',num2str(theta*1000),' mV. Refractory period ',num2str(1000*arp),' ms']);
    xlabel('time (s)')
    ylabel(['Spike threshold (U), number of spikes=',num2str(n_spikes)])
    grid on