import pandas as pd
import numpy as np
from scipy.stats import linregress
from matplotlib import pyplot as plt

## A lot of work to be done here.
## For now, a collection of functions and classes to analyze qPCR data.

#### Intended workflow:
"""
1. Load qPCR data from CSV file into a DataFrame.
2. Identify standards and samples.
3. For each target, fit a standard curve (log concentration vs Cp).
4. For each sample, calculate concentration using the standard curve.
5. For copy numbers calculate ratios (and their standard deviations).
6. Save standard curves as png, save results to a new CSV file.
"""


#### Random functions

# Calculate standard deviation for kan_per_genome using error propagation
# Relative error propagation formula is:
# (Δz/z)² = (Δx/x)² + (Δy/y)² for z = x/y
# Therefore, Δz = z * sqrt((Δx/x)² + (Δy/y)²)

def calculate_stdev_ratio(x, y, z, x_stdev, y_stdev):
    z_stdev = z * np.sqrt((x_stdev / x) ** 2 + (y_stdev / y) ** 2)
    return z_stdev


#### Classes

class Standard:

    def __init__(self, target, df):
        self.target = target
        self.data = df[(df['Target'] == target) & (df['Standard'])]
        self.slope, self.intercept, self.r2 = self.calculate_standard_curve()
        self.best_slope, self.best_intercept, self.best_r2, self.used_data = self.best_standard_curve()


    def calculate_standard_curve(self):
        # Perform linear regression on the log of concentration vs Cp
        log_conc = np.log10(self.data['Concentration'])
        slope, intercept, r2, _, _ = linregress(self.data['Cp'], log_conc)
        return slope, intercept, r2


    def best_standard_curve(self, min_improvement=0.01):
        """
        Try excluding up to 1/3 of points from the start or end.
        Only exclude if R² improves by at least min_improvement.
        Returns: slope, intercept, used_data
        """
        n = len(self.data)
        max_exclude = n // 3
        best_r2 = -np.inf
        best_result = None

        # Prepare data
        sorted_data = self.data.sort_values('Cp')
        log_conc = np.log10(sorted_data['Concentration'].values)
        cp = sorted_data['Cp'].values

        # Try all exclusions from start and end
        for exclude_start in range(0, max_exclude + 1):
            for exclude_end in range(0, max_exclude + 1):
                if exclude_start + exclude_end > max_exclude:
                    continue
                cp_sub = cp[exclude_start:n-exclude_end]
                log_conc_sub = log_conc[exclude_start:n-exclude_end]
                if len(cp_sub) < 3:
                    continue  # Need at least 3 points for regression
                slope, intercept, r_value, _, _ = linregress(cp_sub, log_conc_sub)
                if r_value**2 > best_r2:
                    best_r2 = r_value**2
                    best_result = (slope, intercept, exclude_start, exclude_end, cp_sub, log_conc_sub)

        # Compare with no exclusion
        slope_all, intercept_all, r_value_all, _, _ = linregress(cp, log_conc)
        r2_all = r_value_all**2

        # Only exclude if improvement is substantial
        if best_r2 - r2_all >= min_improvement:
            slope, intercept, exclude_start, exclude_end, cp_sub, log_conc_sub = best_result
            used_data = sorted_data.iloc[exclude_start:n-exclude_end]
            print(f"Excluded {exclude_start} from start and {exclude_end} from end for better fit (R²: {best_r2:.3f} vs {r2_all:.3f})")
            return slope, intercept, best_r2, used_data
        else:
            print("No substantial improvement by excluding points.")
            return slope_all, intercept_all, r2_all, sorted_data
        

    def plot_standard_curve(self, save=True, save_path=None):
        plt.scatter(self.data['Cp'], self.data["Concentration"])
        plt.yscale('log')

        # Plot the fitted standard curve for kan
        cp_range = np.linspace(self.used_data["Cp"].min(), self.used_data["Cp"].max(), 100)
        fitted_conc = 10 ** (self.best_slope * cp_range + self.best_intercept)
        plt.plot(cp_range, fitted_conc, color='red', label=f'Standard Curve ({self.target})')

        plt.xlabel("Cp")
        plt.ylabel("Concentration")
        plt.legend()
        plt.title(f"{self.target} Standard Curve")

        if save:
            plt.savefig(save_path)

        plt.show()


class Sample:

    def __init__(self, name, target, cp: np.array, standard: Standard):
        self.name = name
        self.target = target
        self.cp = cp
        self.standard = standard
        self.concentration = self.calculate_concentration()
        self.mean_conc = np.mean(self.concentration)
        self.std_conc = np.std(self.concentration)

    def calculate_concentration(self):
        slope = self.standard.best_slope
        intercept = self.standard.best_intercept
        conc = 10 ** (slope * self.cp + intercept)
        return conc


working_dir = '/Users/nataschaspahr/Downloads/qPCR_Data_Natascha/qPCR_analysis/'
qPCRfiles = [
    "ALE1b_firstTransfer_pyruvate_5-28-25_qPCR.csv",
    "ALE1b_firstTransfer_succinate_6-24-25_qPCR.csv",
    "ALE1b_lastTransfer_pyruvate_5-28-25_qPCR.csv",
    "ALE1b_lastTransfer_succinate_7-10-25_qPCR.csv",
    "ALE1b_overnightCultures_5-21-25_qPCR.csv"
]

for file in qPCRfiles:

    print(file)
    # Load the data into a DataFrame
    df = pd.read_csv(working_dir + file)[['Name', 'Standard', 'Target', 'Concentration', 'Cp']].convert_dtypes()
    
    # Ensure that all standards have concentrations
    try:
        assert not df.loc[df['Standard']]['Concentration'].isnull().any(), "Some standards are missing concentrations!"
    except AssertionError as e:
        print(e)

    # Determine number of targets
    targets = df['Target'].dropna().unique()
    print(f"Targets: {targets}")

    # Ensure that each target has its own standard

    def check_targets_have_standards(df):
        """
        Checks if each target in the DataFrame has at least one associated standard.
        Prints a warning for any target without a standard.
        """
        targets = df['Target'].dropna().unique()
        missing = []
        for target in targets:
            # Check if any row for this target is marked as a standard
            has_standard = ((df['Target'] == target) & (df['Standard'])).any()
            if not has_standard:
                missing.append(target)
        if missing:
            print(f"Warning: The following targets have no associated standard: {missing}")
        else:
            print("All targets have at least one associated standard.")

    check_targets_have_standards(df)

    kan_stan = Standard('kan', df)
    kan_stan.plot_standard_curve(save=True, save_path=working_dir+file.replace(".csv", "_kan_standard.png"))
                                 
    samples = df.loc[~df['Standard'] & ~df['Name'].isnull()]

    results = []

    for name, group in samples.groupby('Name'):
        kan_sample = Sample(name, 'kan', group.loc[group['Target']=='kan']['Cp'].to_numpy(), kan_stan)
        rpoA_sample = Sample(name, 'rpoA', group.loc[group['Target']=='rpoA']['Cp'].to_numpy(), kan_stan)
        kan_rpoA_ratio = kan_sample.mean_conc / rpoA_sample.mean_conc
        ratio_std = calculate_stdev_ratio(
            kan_sample.mean_conc, rpoA_sample.mean_conc, kan_rpoA_ratio,
            kan_sample.std_conc, rpoA_sample.std_conc
            )
        results.append({
            'Sample': name,
            'Kan_Mean_Conc': kan_sample.mean_conc,
            'Kan_Std_Conc': kan_sample.std_conc,
            'rpoA_Mean_Conc': rpoA_sample.mean_conc,
            'rpoA_Std_Conc': rpoA_sample.std_conc,
            'kan_rpoA_Ratio': kan_rpoA_ratio,
            'kan_rpoA_Ratio_Std': ratio_std
        })

    out_frame = pd.DataFrame(results)
    out_frame.to_csv(working_dir + file.replace(".csv", '_qPCR_results.csv'), index=False)