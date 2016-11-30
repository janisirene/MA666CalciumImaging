%% Autoregressive Model for Calcium Dynamics
% Explain what an autoregressive model is and how this applies

%% Generate spikes from exponential distribution
rate = 1000; % spike rate (Hz)
dur = .01; % duration of trial (s)
dt = 0.0001; % time step (s)
T = (0:dt:dur);

spikes = zeros(dur/dt+1,1);

% Generating spikes from a exponential distribution
for t=1:length(T)
    if (rate*dt)>=rand
        spikes(t) = 1;
    end
end

%% Applying the AR model to Action Potentials
% Before implementing the calculation, consider what each of the parameters
% in the model does. Look again at the AR model.
%%
% $c(t) = \sum_{j=1}^p\gamma_jc(t-j)+s(t)$ 
%%
% Consider the following for $p=1$:
%%
% * What does the parameter $\gamma_1$ represent? 
% * How do changes to $\gamma_1$ affect the calcium signal? 
% * What are reasonable values for this parameter?
%%
% Look at what happens when we vary $\gamma_1$:
p = 1;
Gamma = (.5:.05:.99);
calcium = zeros(length(spikes),numel(Gamma));
for g = 1:numel(Gamma)
    gamma = Gamma(g);
    for t = 2:length(spikes)
        calcium(t,g) = gamma(p)*calcium(t-p,g)+spikes(t);
    end
    figure(1); plot((0:dt:dur),calcium(:,g)+2*g); hold on
end
title('Calcium concentration for varying \gamma')
hold off;
legend(num2str(Gamma'))
%%
% Can you tell what's happening?
% What if we overlay the plots of the calcium signal with the highest and
% lowest values of $\gamma$.
figure(2); plot((0:dt:dur),calcium(:,[1,g]))
legend(num2str(Gamma([1,g])'))
title('Comparison of traces with high and low \gamma')
%%
% Let's do one more comparison: are the peaks in the same place?

[~,peaklocs] = arrayfun(@(i) findpeaks(calcium(:,i)),[1;g],...
    'uniformoutput',false);
disp(cat(2,peaklocs{:}))

%%
% Great! Looks like yes! Now, let's look at what happens when $p=2$.
p = 2;
gamma1 = 0.5;
Gamma = (0:.05:.49);
calcium = zeros(length(spikes),numel(Gamma));
for g = 1:numel(Gamma)
    gamma = [gamma1; Gamma(g)];
    for t = 3:length(spikes)
        calcium(t,g) = gamma(1)*calcium(t-1,g)+gamma(2)*calcium(t-2,g)+spikes(t);
%         calcium(t,g) = sum(gamma.*calcium(t-(1:p),g))+spikes(t);
%         calcium(t,g) = gamma*calcium(t-p,g)+spikes(t);
    end
    figure(3); plot((0:dt:dur),calcium(:,g)+2*g); hold on
end
title('Calcium concentration for \gamma_1=0.5 with varying \gamma_2')
hold off;
legend(num2str(Gamma'))

%%
% Look at the values of |Gamma|. Why do you think those values were chosen?
% What happens if you let |Gamma| go to 0.7 for example? Change the
% parameters and see what happens.