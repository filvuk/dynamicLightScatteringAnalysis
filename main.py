from pathlib import Path
from collections import defaultdict


def load_data(folder: Path):
    runs = (run for run in folder.iterdir() if run.suffix.lower() == ".asc")
    experiment_name = folder.stem
    # dict containing data of runs. Auto-sorts repetitions under same key
    runs_data = defaultdict(lambda: defaultdict(list))
    for run in runs:
        # get run number from file name
        run_number = int(run.stem[-4])
        with open(run, "r") as content:
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
            angle_rounding_prec = 1
            assert "Angle" in lines[18].split(":")[0]
            properties["angle"] = round(float(lines[18].split(":")[1].strip()), angle_rounding_prec)

            # create dictionary hash to sort by repetitions
            dict_hash = hash(f"{properties["temperature"]}{properties["viscosity"]}{properties["refractive_index"]}{properties["angle"]}")
            runs_data[dict_hash].update(properties)
            assert "Correlation" in lines[25]
            for corr_line in lines[26:]:
                corr_line.
            runs_data[dict_hash]["runs"].append()






def plot_radius_distribution():
    data = load_data()
    result = data.calculate()
    plot(result)
