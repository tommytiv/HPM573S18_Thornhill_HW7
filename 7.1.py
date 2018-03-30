from enum import Enum
import numpy as np


class HealthStat(Enum):
    """ health status of patients  """
    ALIVE = 1
    DEAD = 0


class Patient(object):
    def __init__(self, id, mortality_prob):
        """ initiates a patient
        :param id: ID of the patient
        :param mortality_prob: probability of death during a time-step (must be in [0,1])
        """
        self._id = id
        self._rnd = np.random       # random number generator for this patient
        self._rnd.seed(self._id)    # specifying the seed of random number generator for this patient

        self._mortalityProb = mortality_prob
        self._healthState = HealthStat.ALIVE  # assuming all patients are alive at the beginning
        self._survivalTime = 0

    def simulate(self, n_time_steps):
        """ simulate the patient over the specified simulation length """

        t = 0  # simulation current time

        # while the patient is alive and simulation length is not yet reached
        while self._healthState == HealthStat.ALIVE and t < n_time_steps:
            # determine if the patient will die during this time-step
            if self._rnd.sample() < self._mortalityProb:
                self._healthState = HealthStat.DEAD
                self._survivalTime = t + 1  # assuming deaths occurs at the end of this period

            # increment time
            t += 1

    def get_survival_time(self):
        """ returns the patient survival time """
        # return survival time only if the patient has died
        if self._healthState == HealthStat.DEAD:
            return self._survivalTime
        else:
            return None


class Cohort:
    def __init__(self, id, pop_size, mortality_prob):
        """ create a cohort of patients
        :param id: cohort ID
        :param pop_size: population size of this cohort
        :param mortality_prob: probability of death for each patient in this cohort over a time-step (must be in [0,1])
        """
        self._patients = []      # list of patients
        self._survivalTimes = []  # list to store survival time of each patient

        # populate the cohort
        for i in range(pop_size):
            # create a new patient (use id * pop_size + n as patient id)
            patient = Patient(id * pop_size + i, mortality_prob)
            # add the patient to the cohort
            self._patients.append(patient)

    def simulate(self, n_time_steps):
        """ simulate the cohort of patients over the specified number of time-steps
        :param n_time_steps: number of time steps to simulate the cohort
        """
        # simulate all patients
        for patient in self._patients:
            # simulate
            patient.simulate(n_time_steps)
            # record survival time
            value = patient.get_survival_time()
            if not (value is None):
                self._survivalTimes.append(value)

    def get_ave_survival_time(self):
        """ returns the average survival time of patients in this cohort """
        return sum(self._survivalTimes)/len(self._survivalTimes)

    def get_over_time(self):
        count_over = 0
        for x in self._survivalTimes:
            if x > YEAR:
                count_over += 1
        return count_over


MORTALITY_PROB = 0.1    # annual probability of mortality
TIME_STEPS = 100        # simulation length
SIM_POP_SIZE = 573     # population size of the simulated cohort
ALPHA = 0.05            # significance level
YEAR = 5


# create a cohort of patients
myCohort = Cohort(id=1, pop_size=SIM_POP_SIZE, mortality_prob=MORTALITY_PROB)

# simulate the cohort
myCohort.simulate(TIME_STEPS)

# print the patient survival time
print('Average survival time (years):', myCohort.get_ave_survival_time())
#print('95% CI of average survival time (years)', myCohort.get_CI_survival_time(ALPHA))

print ("Patients who survived over 5 years: ", myCohort.get_over_time())
print ("Percentage of patients who survived over 5 years: ", myCohort.get_over_time()/SIM_POP_SIZE)
