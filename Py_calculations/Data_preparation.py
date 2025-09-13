

import pandas as pd
import numpy as np
import pickle
import random

# read data
df = pd.read_csv("wheeler.csv", sep = '\t')

# fix column names
better_names = [str.strip(c) for c in df.columns]
df.columns = better_names

# choose relevant catalysts
good_cats = ["1a", "1b", "2", "3a", "3b", "4a", "4b", "4c", "4d",
             "4e", "5a", "5b", "5c", "6"]
df = df[df['Catalyst'].isin(good_cats)]


# CONSTANTS
HARTREE_TO_KCAL_PER_MOLE = 627.5095


# Define the number of TS energies to keep and the large number to fill gaps
N = 10
Nbig = 100.0

# Get the list of methods
methods = df.columns[3:]

# Initialize dictionaries to store the matrices
R_matrices = {}
S_matrices = {}

# Iterate over each method
for method in methods:
    # Filter the DataFrame for the current method
    method_df = df[['Catalyst', 'Product', method]]

    # Get the list of unique catalysts
    catalysts = method_df['Catalyst'].unique()

    # Initialize the R and S matrices with Nbig values
    R_matrix = np.full((len(catalysts), N), Nbig)
    S_matrix = np.full((len(catalysts), N), Nbig)

    # Iterate over each catalyst
    for i, catalyst in enumerate(catalysts):
        # Filter the DataFrame for the current catalyst
        catalyst_df = method_df[method_df['Catalyst'] == catalyst]

        # Get the R and S energies for the current catalyst
        R_energies = catalyst_df[catalyst_df['Product'] == 'R'][method].values
        S_energies = catalyst_df[catalyst_df['Product'] == 'S'][method].values

        # Sort the energies in ascending order and delete nans
        R_energies.sort()
        S_energies.sort()

        R_energies = R_energies[~np.isnan(R_energies)]
        S_energies = S_energies[~np.isnan(S_energies)]

        # Substract minimum value and transform from Hartree to kcal/mol
        min_E = np.min(np.concatenate([R_energies,
                                       S_energies]))
        R_energies = ( R_energies - min_E ) * HARTREE_TO_KCAL_PER_MOLE
        S_energies = ( S_energies - min_E ) * HARTREE_TO_KCAL_PER_MOLE

        # Fill the R and S matrices with the sorted energies (up to N values)
        R_matrix[i, :len(R_energies[:N])] = R_energies[:N]
        S_matrix[i, :len(S_energies[:N])] = S_energies[:N]

    # Store the matrices in the dictionaries
    R_matrices[method] = R_matrix
    S_matrices[method] = S_matrix


exp_ee  = [74, 72, 88, 52, 81, 65, -48,
           46, -53, 44, 80, 61, -49, 84]

temp = [195, 195, 195, 296, 193, 233, 233,
        233, 195, 195, 195, 195, 195, 228]

cats = ['1a', '1b', '2', '3a', '3b', '4a', '4b',
        '4c', '4d', '4e', '5a', '5b', '5c', '6']


with open("wheeler.pcl", 'wb') as f:
    pickle.dump({'Methods': methods, 'Rs': R_matrices,
                 'Ss': S_matrices, 'Ts': temp, 'EEs': exp_ee }, f)


# In order to delete random lowest TS (R or S) add these lines before substracting minimum value and transforming from Hartree to kcal/mol:
# # Find the two minimum energies: one from R and one from S
#     min_R = np.min(R_energies)
#     min_S = np.min(S_energies)
# # Randomly decide whether to remove the minimum from R or S
#     if random.choice([True, False]):
#         R_energies = R_energies[R_energies != min_R]
#     else:
#         S_energies = S_energies[S_energies != min_S]






