import scipy.stats as stat


SIM_POP_SIZE = 1000       # population size of simulated cohorts
TIME_STEPS = 100        # length of simulation
ALPHA = 0.05             # significance level for calculating confidence intervals
NUM_SIM_COHORTS = 1000   # number of simulated cohorts used to calculate prediction intervals
PROBABILITY = 0.50

# details of a clinical study estimating the mean survival time
OBS_N = 573        # number of patients involved in the study
OBS_MEAN = 400/573    # estimated mean survival time
#OBS_HL =       # half-length
OBS_ALPHA = 0.05   # significance level
# the standard deviation of the mean survival time reported in the clinical study
# assumes that the reported confidence interval in this study is a t-confidence interval
#OBS_STDEV =  #OBS_HL / stat.t.ppf(1 - OBS_ALPHA / 2, OBS_N-1)

# how to sample the posterior distribution of mortality probability
# minimum, maximum and the number of samples for the mortality probablity
POST_L, POST_U, POST_N = 0.45, 0.75, 1000
