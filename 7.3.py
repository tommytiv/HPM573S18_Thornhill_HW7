import scipy.stats as stat

# Parameters for a binomial probability mass function
x = 400     # x = total number of “successes”
n = 573     # n = number of trials
p = 0.5    # p = probability of a success on an individual trial

Probability = stat.binom._pmf(x,n,p)
SD = stat.binom.std(n, p, loc=0)

print ("The likelihood that a clinical study reports 400 of 573 "
      "participants survived at the end of the 5-year study period "
       "if 50% of the patients in our simulated cohort survived beyond 5 years is ")
print (Probability)

