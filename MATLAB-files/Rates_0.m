function Rate = Rates_0(Spikes,window_size,bin_size)
%Transforming spike trains (binary series) to single-trial time-dependent firing rates
%by a smoothign rectangular time window.
%Emili B-B, Bournemouth University. 
%Neuronal Analysis Masters Course.
%Suboptimal code, desgined for education/portability to other languages.
%Optional inputs: 
%           Spikes=Matrix, each column is a neuron, each row is an
%                       observation. Last column contains the trial number.
%           window_size=Rectangular window size in seconds, see definition of time-dependent firing rate (or just rate).
%           bin_size=max_resolution of each matrix row, in seconds (must match the
%               criterion used in input #1, matrix 'Spikes').
%Output: 
%          Rate=Matrix of length number of rows in the input matrix
%          'Spikes' and the same number of columns as
%          the input matrix 'Spikes'. Each column is the
%          time-dependent firing-rate in spikes/second. 
%          The last column indicates the trial number.
%Last modification: 19/10/2022. 
%% PARAMETER SETUP

%1-Parameters to modify
if nargin<1    
    load('DataSetQ3.txt'); 
    Spikes=DataSetQ3;
end
if nargin<2
   window_size=0.2;%In sec.This is the period of study/sliding time window we will use to sum spikes
end   
%2-Fixed parameter not to modify, since it depends on the 'spikes' rows
%definition.
if nargin<3
    bin_size=0.01;%In sec. (10 ms hard-coded)
 end
     
%%   RATES ESTIMATION

%1-Basic trial-to-trial rate estimation
trial_num=Spikes(:,end);
spikes=Spikes(:,1:end-1); %Trial number removed
[n,m]=size(spikes);
times=[0:bin_size:(n-1)*bin_size];
window_steps=round(window_size/bin_size);%Number of steps in a time window 
Rate=zeros(n,m);
%Two nested loops: over neurons and for each spike train
for j=1:m %Neurons index. A loop is used just for clarity/portability to other languages even thou is much slower than a more compact implementation.
    %this is the same as "for j in range(m)" in python
    Smoothed=[];
    for i=1:(n-window_steps+1) %Suggest to check why the "+1": substitute this latest value of the loop in the line below
        Spikes_one_neuron=spikes(i:i+window_steps-1,j);%Selecting a column i.e., a neuron
        Rate_in_this_window=sum(Spikes_one_neuron).*(1/window_size); %Just the definition of 
        %Firing-rate in spikes/s. Note the "-1" in the slice of the matrix, since we start the count on "+1"
        Smoothed=[Smoothed,Rate_in_this_window];%Simply stacking rate matrices
    end
    Smoothed=[Smoothed,repmat(Smoothed(end),1,window_steps-1)]';%Complete the series by repeating the last rate value until the "Smoothed" vector is of length "n" 
    Rate(:,j)=Smoothed;
end

Rate=[Rate,trial_num]; %Adding back the trial number

 
 %% PLOTTING
 figure 
 for j=1:m %Neurons
    subplot (m,2,2*j-1)
    plot(times,Spikes(:,j))
    title(['Neuron ',num2str(j)'])
    ylabel('Spikes')
end
xlabel('Time (s)')

for j=1:m %Neurons
    subplot (m,2,2*j)
    plot(times,Rate(:,j),'r')
    title(['Neuron ',num2str(j),' rate'])
    ylabel('Spikes/sec')
end
xlabel('Time (s)')

end
