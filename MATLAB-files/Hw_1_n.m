function Rate = Rates_0(Spikes,window_size,bin_size)

if nargin<1
    load('DataSetQ3.txt');
    Spikes=DataSetQ3;
end
if nargin<2
    window_size=0.2;
end
if nargin<3
    bin_size=0.01;
end

trial_num=Spikes(:,end);
spikes=Spikes(:,1:end-1);
[n,m]=size(spikes);
times=[0:bin_size:(n-1)*bin_size];
window_steps=round(window_size/bin_size);
Rate=zeros(n,m);

for j=1:m
    Smoothed=[];
    for i=1:(n-window_steps+1)
        Spikes_one_neuron=spikes(i:i+window_steps-1,j);
        Rate_in_this_window=sum(Spikes_one_neuron).*(1/window_size);
        Smoothed=[Smoothed,Rate_in_this_window];
    end
    Smoothed=[Smoothed,repmat(Smoothed(end),1,window_steps-1)]';
    Rate(:,j)=Smoothed;
end

Rate=[Rate,trial_num];

figure
for j=1:m
    subplot(m,2,2*j-1)
    plot(times,Spikes(:,j))
    title(['Neuron',num2str(j)'])
    ylabel('Spikes')
end
xlabel('Time(s)')

for j=1:m
    subplot(m,2,2*j)
    plot(times,Rate(:,j),'r')
    title(['Neuron',num2str(j),'rate'])
    ylabel('Spikes/sec')
end
xlabel('Times(s)')

end

