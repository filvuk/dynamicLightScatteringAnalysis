from pathlib import Path
from collections import defaultdict
import helpers
import numpy as np
from matplotlib import pyplot as plt
from scipy.signal import argrelmax


def load_data(folder: Path, validation: bool):
    runs = (run for run in folder.iterdir() if run.suffix.lower() == ".asc")
    # dict containing data of runs. Auto-sorts repetitions under same key
    runs_data = defaultdict(lambda: defaultdict(list))
    runs_data["experiment_name"] = folder.stem
    for run in runs:
        # get run number from file name
        run_number = int(run.stem[-4])
        with open(run, "r", encoding="iso-8859-1") as content:
            lines = content.readlines()
            properties = {}

            temp_rouding_prec = 1
            assert "Temperature" in lines[14].split(":")[0]
            properties["temperature"] = round(float(lines[14].split(":")[1].strip()), temp_rouding_prec)
            visc_rounding_prec = 2
            assert "Viscosity" in lines[15].split(":")[0]
            properties["viscosity"] = round(float(lines[15].split(":")[1].strip()), visc_rounding_prec)
            refr_rounding_prec = 2
            assert "Refractive Index" in lines[16].split(":")[0]
            properties["refractive_index"] = round(float(lines[16].split(":")[1].strip()), refr_rounding_prec)
            lambda_rounding_prec = 0
            assert "Wavelength" in lines[17].split(":")[0]
            properties["wavelength"] = round(float(lines[17].split(":")[1].strip()), lambda_rounding_prec)
            angle_rounding_prec = 1
            assert "Angle" in lines[18].split(":")[0]
            properties["angle"] = round(float(lines[18].split(":")[1].strip()), angle_rounding_prec)

            # create dictionary hash to sort by repetitions
            dict_hash = hash(f"{properties["temperature"]}{properties["viscosity"]}{properties["refractive_index"]}{properties["wavelength"]}{properties["angle"]}")
            runs_data[dict_hash].update(properties)
            runs_data[dict_hash]["run_numbers"].append(run_number)

            assert "Correlation" in lines[25]
            data = defaultdict(list)
            for corr_line in lines[26:]:
                if not corr_line.strip():
                    break
                values = corr_line.split()
                if not validation:
                    data["times"].append(float(values[0]) / 1e3)
                    data["ys"].append(float(values[1]) + 1)
                else:
                    data["times"].append(float(values[0]) / 1e6)
                    data["ys"].append(float(values[1]))
            runs_data[dict_hash]["runs"].append(data)
    return runs_data

def calculate(data):
    for experiment in data.values():
        if isinstance(experiment, str):
            continue
        refractive_index = complex(experiment["refractive_index"]).real
        wavelength = experiment["wavelength"]
        angle = experiment["angle"] / 180 * np.pi
        temperature = experiment["temperature"]
        viscosity = experiment["viscosity"] * 1e-3 # conversion to pascal seconds from centipoise
        
        lowHr = 0.09
        highHr = 1e6
        n = 200
        q = helpers.get_q(wavelength, refractive_index, angle)
        sUpLimitHigh = helpers.s_inverse_decay_rate(
            helpers.diffusion_from_hydrodynamic_radius(
                highHr/1e9,
                temperature,
                viscosity
            ),
            q
        )

        sUpLimitLow = helpers.s_inverse_decay_rate(
            helpers.diffusion_from_hydrodynamic_radius(
                lowHr/1e9,
                temperature,
                viscosity
            ),
            q
        )

        # Sequence in linear space! 10.0**start to 10**stop
        s_space = np.logspace(np.log10(sUpLimitLow),np.log10(sUpLimitHigh), n) 

        ds = helpers.diffusion_from_inverse_decay_rate(s_space,q)
        hrs = helpers.hydrodynamic_radius(ds, temperature, viscosity)*1e9  # In nanometers
        for run in experiment["runs"]:
            times = np.array(run["times"])
            ys = np.array(run["ys"])
            betaGuess = helpers.get_beta_prior_single(ys, times) 
            g1 = np.array(helpers.g1_from_g2(ys, betaGuess))

            alpha = 0.1
            timeLimit = 1e8
            selectedTimes = times < (timeLimit / 1e6)
            # Return the fitted contributions and residuals of the first order autocorrelation function
            contributionsGuess, residualsG1, _ = helpers.get_contributios_prior_single(g1[selectedTimes], times[selectedTimes], s_space, betaGuess, alpha)

            curvesResidualNorm, curvesPenaltyNorm = [],[]
            alphaVec = (5**np.arange(-6,2,0.1,dtype=float))**2
            # Iterate over the vector with different values of alpha
            for alpha in alphaVec:
                _ , residualNorm, penaltyNorm = helpers.get_contributios_prior_single(
                    g1[selectedTimes],times[selectedTimes],
                    s_space,betaGuess,alpha) 
                curvesResidualNorm.append(residualNorm) # List (one element per alpha)
                curvesPenaltyNorm.append(penaltyNorm)   # List (one element per alpha)
            curvesResidualNorm = np.array([i[0] for i in curvesResidualNorm]) # One row per alpha
            curvesPenaltyNorm = np.array([i[0] for i in curvesPenaltyNorm])  # One row per alpha

            alphaOptIdx = helpers.find_Lcurve_corner_single(curvesResidualNorm,curvesPenaltyNorm)
                
            optimalAlpha = alphaVec[alphaOptIdx]
            contributionsGuess, residualsG1, _ = helpers.get_contributios_prior_single(g1[selectedTimes], times[selectedTimes], s_space, betaGuess, optimalAlpha)
            curvesResidualNorm, curvesPenaltyNorm = [],[]
            # Iterate over the vector with different values of alpha
            for alpha in alphaVec:
                _ , residualNorm, penaltyNorm = helpers.get_contributios_prior_single(
                    g1[selectedTimes],times[selectedTimes],
                    s_space,betaGuess,alpha) 
                curvesResidualNorm.append(residualNorm) # List (one element per alpha)
                curvesPenaltyNorm.append(penaltyNorm)   # List (one element per alpha)
            curvesResidualNorm = np.array(curvesResidualNorm) # One row per alpha
            curvesPenaltyNorm  = np.array(curvesPenaltyNorm)  # One row per alpha
            run["hdr_dist"] = (hrs, contributionsGuess[0])
    

def plot_radius_distribution(folder: Path, validation = False):
    data = load_data(folder, validation)
    calculate(data)
    experiment_name = ""
    for experiment in data.values():
        if isinstance(experiment, str):
            experiment_name = experiment
            continue
        wavelength = experiment["wavelength"]
        angle = experiment["angle"]
        temperature = experiment["temperature"]
        for run in experiment["runs"]: 
            plt.xscale("log")
            plt.plot(run["hdr_dist"][0],run["hdr_dist"][1])
            plt.xlabel("Hydrodynamic radius (nm)")
            plt.ylabel("Relative contribution")
            plt.title(fr"{experiment_name}: HDR distribution at $\lambda = ${wavelength} nm, {angle} °, {temperature} K")
            plt.show()
            x_max = argrelmax(run["hdr_dist"][1], order = 5)
            print(f"Particle size at maximum is {np.round(run["hdr_dist"][0][x_max], 1)} nm")
