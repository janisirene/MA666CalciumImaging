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
% * What does the parameter $\gamma_1$ represent? 
% * How do changes to $\gamma_1$ affect the calcium signal? 
% * What are reasonable values for this parameter?
%%
% Look at what happens when we vary $\gamma_1$:
P = 1;
Gamma = (.5:.05:.99);
calcium = zeros(length(spikes),numel(Gamma));
for g = 1:numel(Gamma)
    gamma = Gamma(g);
    for t = 2:length(spikes)
        for p = 1:P
            calcium(t,g) = calcium(t) + gamma(p)*calcium(t-p)+spikes(t);
        end
    end
    plot((0:dt:dur),calcium(:,g)+2*g); hold on
end
hold off;
legend(num2str(Gamma'))
%%
% Can you tell what's happening?
% What if we overlay the plots of the calcium signal with the highest and
% lowest values of $\gamma$.
plot((0:dt:dur),calcium(:,[1,g]))
legend(num2str(Gamma([1,g])'))