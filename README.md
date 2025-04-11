[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/ale-tom/PD_HDDM/HEAD)
# PD_HDDM 

## Repository description
The code in this repository reproduces most analyses and figures from Tomassini et al.,2019 **'Learning from the past and expecting the future in Parkinsonism: Dopaminergic influence on predictions about the timing of future events'** published in Neuropsychologia: https://doi.org/10.1016/j.neuropsychologia.2019.02.003.
Model fitting is implemented on iPython using the Hierarchical Drift-Diffusion Modelling toolbox (https://github.com/hddm-devs/hddm). Each model tested in the paper is implemented on a separate Jupyter notebook. Please note that model fitting will require a long time (approx 15 hours on a Mac Pro with 4 CPUs). Estimated model parameters (including DICs) are reported in Stats_Model#.txt files.

## Rationale and findings
Parkinson’s disease slows down movement and reaction times, mainly because the brain doesn’t have enough dopamine. Dopamine helps us use past experiences to predict what will happen next. Our research shows that when dopamine is low, people with Parkinson’s have trouble predicting when to move, making their timing less accurate and more uncertain. Importantly, we found that medication that restores dopamine helps fix this problem making their movements faster and more timely.

<br>
<br>
<br>
<br>

![1-s2 0-S0028393218306171-gr2_lrg](https://github.com/user-attachments/assets/46d25814-f68a-4b12-b22e-dc87e6ed6c19)
*Behavioural results (a) Modulation of RTs by temporal predictability. Average RTs are plotted against foreperiod variability (bottom abscissa) and segregated by foreperiod duration (upper abscissa) in controls (green), PD-on (red) and PD-off (blue) patients. (b) Poisson regression. The relationship between RTs and temporal predictability were characterized in terms of intercept (β0) and slope (β1) of a linear function fitted to the RTs of each participant. * p < 0.05; Error bars represent ± 1 SEM. (c) Foreperiod analysis. RTs were averaged across subjects and conditions and plotted as a function of foreperiod duration. RTs were normalized and smoothed only for illustrative purposes. The inset barplot shows the slopes of the linear regression between RTs and foreperiod duration. Negative slopes indicate a normal foreperiod effect. No significant difference between groups was detected.*

<br>
<br>
<br>


![1-s2 0-S0028393218306171-gr3_lrg](https://github.com/user-attachments/assets/8d348fe7-8f91-4378-8596-7c641560c065)
*Model comparison and hierarchical drift-diffusion modelling of variable foreperiod task (a) Model comparison. The deviance information criterion (DIC) between the best fitting model (Model 1 - free parameter: boundary) and the two alternative models (Model 2 - drift-rate; Model 3 – boundary and drift rate) is shown. (b) Influence of predictability on evidence accumulation. Intercept and slope of a linear function were fitted to the estimated boundary of each participant. (c) Model fit. Estimated boundary values were strongly correlated to empirical RT averaged across trials within each predictability level (Spearman correlation); * p < 0.05; Error bars represent ± 1 SEM.*

