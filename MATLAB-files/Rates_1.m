function Rate = Rates_1(Spikes,window_size,bin_size)
%Transforming spike trains (binary series) to single-trial time-dependent
%firing rates via a rectangular window and comparing them using a correlation coefficent.
%Emili B-B, Bournemouth University. 
%Neuronal Analysis Masters Course. 
%Suboptimal code, desgined for pedagogical pruposes/portability to other languages.
%Optional inputs: 
%           Spikes=Matrix, each column is a neruon, each row is an observation.
%                        Last column contains the trial number.
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
figure
%1-Parameters to modify
if nargin<1    
    load('DataSetQ3.txt'); 
    Spikes=DataSetQ3;
end
if nargin<2
   window_size=0.5;%In sec. This is the period of study/sliding time window we will use to sum spikes
end   

%2-Fixed parameter not to modify, since it depends on the 'spikes' rows
%definition (10 ms hard-coded).
if nargin<3
    bin_size=0.01;%In sec.
end
 

     
%%   RATES ESTIMATION

%1-Basic trial-to-trial rate estimation
trial_num=Spikes(:,end);
spikes=Spikes(:,1:end-1); %Trial number removed
[n,m]=size(spikes);
times=[0:bin_size:(n-1)*bin_size];
window_steps=round(window_size/bin_size);%number of steps in a time window 
Rate=zeros(n,m);
%Two nested loops: over neurons and for each spike train
for j=1:m %Neurons index. A loop is used just for clarity/portability to other languages even thou is much slower than a more compact implementation.
    %this is the same as "for j in range(m)" in python
    Smoothed=[];
    for i=1:(n-window_steps+1)  %Suggest to check why the "+1": substitute this latest value of the loop in the line below   
        Spikes_one_neuron=spikes(i:i+window_steps-1,j);%Selecting a column i.e., a neuron
        Rate_in_this_window=sum(Spikes_one_neuron).*(1/window_size); %Just the definition of 
        %Firing-rate in spikes/s. Note the "-1" in the slice of the matrix, since we start the count on "+1"
        Smoothed=[Smoothed,Rate_in_this_window];%Simply stacking rate matrices   
        % Smoothed=[Smoothed,( sum(spikes((i:i+window_steps-1),j)).*(1/window_size) )]; %same but more compact
    end
    Smoothed=[Smoothed,repmat(Smoothed(end),1,window_steps-1)]';%Complete the series by repeating the last rate value until the "Smoothed" vector is of length "n" 
    Rate(:,j)=Smoothed;
end

Rate=[Rate,trial_num]; %Adding back the trial number

 
 %% PLOTTING
 
 for j=1:m %Neurons
    subplot (m,2,2*j-1)
    plot(times,Spikes(:,j))
    title(['Neuron ',num2str(j)'])
    ylabel('Spikes'),
end
xlabel('Time (s)')

for j=1:m %Neurons
    subplot (m,2,2*j)
    plot(times,Rate(:,j),'r')
    title(['Neuron ',num2str(j),' rate'])
    ylabel('Spikes/sec')
end
xlabel('Time (s)')



%% ADDITION to Rates_0: COMPARISON 
%Correlation coeffients between spike trains and rates
disp('********************************')
disp('CORRELATION COEFFICIENTS RATE vs SPIKES (p value). Significance: p<0.05');  
disp('********************************')
for j=1:m
      for i=1:max(trial_num)
         trial_index=(trial_num(:,end)==i);%This is a logical array (1s or 0s). It is a trick which will help us to select only data for trial "i" below
        current_trial_rates=Rate(trial_index,j);
        current_trial_spikes=Spikes(trial_index,j);
        %Check the function below in matlab -simply type ">help corrcoeff"
        %or on the top right bar in the IDE
        [corr_value,p]=corrcoef(current_trial_rates,current_trial_spikes);corr_value=corr_value(1,2);p=p(1,2);
        %Only the off-diagonal elements of the two outputs matter (the
        %diagonal elements are self-correlation, =1 always). "p" is the
        %statistical significance, how reliable is the corresponding value. p<0.005
        %usually means that the value is trustable with a probability of
        %99.5%
        disp(['Trial ',num2str(i),' Neuron ',num2str(j),' correlation = ',num2str(corr_value),' (',num2str(p),')']);  
      end
      current_rates=Rate(:,j);
      current_spikes=Spikes(:,j);
      [corr_value,p]=corrcoef(current_rates,current_spikes);corr_value=corr_value(1,2);p=p(1,2);
      disp('********************************')
      disp(['All trials Neuron ',num2str(j),' correlation = ',num2str(corr_value),' (',num2str(p),')']);  
      disp('********************************')
end

end