# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 17:07:45 2024

@author: enf21

This code is associated with, "Sloppy feeding by predators destabilises early animal food webs in the Ediacaran".
To run this code, change the "files" variable to point to the processed data outputs to be compared.
"""

import seaborn
import matplotlib.pyplot as plt
import scipy.stats as sci
import numpy as np
import statistics

pixel_count = 60

plot = True
lock_z_axis = True
plotaxis = 0
plot_grid = 0 # 4, for effect of zooplankton in a predator-filled ecosystem, 1 for zoo into a no-predator ecosystem

summary_one_cube = True
summary_which_cube = 0

output_detritus_curve = True
outputfile = "./run.txt"

files = [".processed_NoPredators.txt",".processed_Predators.txt"]
names = ["NoPredNoZoo", "PredNoZoo"]

Variables = ["DOC biomass", "Autotrophic plankton biomass", "Heterotrophic plankton biomass"]

Variable1 = Variables[0]
Variable2 = Variables[1]
Variable3 = Variables[2]

DOC_biomasses = []
index = 0
while index <= 60:
    DOC_biomasses.append(0)
    index = index + 1
DOC_biomasses[0] = np.log(2.4)
DOC_biomasses[60] = np.log(24000)
index = 0
while index <= 60:
    DOC_biomasses[index] = DOC_biomasses[0] + (index / 60) * (DOC_biomasses[60] - DOC_biomasses[0])
    DOC_biomasses[index] = np.exp(DOC_biomasses[index])
    index = index + 1

DOC_intermediate_biomasses = []
index = 0
while index < 60:
    DOC_intermediate_biomasses.append((DOC_biomasses[index] + DOC_biomasses[index + 1]) / 2)
    index = index + 1

grids = []
index = 0
while index < len(names):
    grids.append([])
    index = index + 1

index0 = 0
while index0 < len(names):
    index1 = 0
    while index1 < pixel_count:
        grids[index0].append([])
        index1 = index1 + 1
    index0 = index0 + 1

index0 = 0
while index0 < len(names):
    index1 = 0
    while index1 < pixel_count:
        index2 = 0
        while index2 < pixel_count:
            grids[index0][index1].append([])
            index2 = index2 + 1
        index1 = index1 + 1
    index0 = index0 + 1

index0 = 0
while index0 < len(names):
    index1 = 0
    while index1 < pixel_count:
        index2 = 0
        while index2 < pixel_count:
            index3 = 0
            while index3 < pixel_count:
                grids[index0][index1][index2].append(0)
                index3 = index3 + 1
            index2 = index2 + 1
        index1 = index1 + 1
    index0 = index0 + 1

n = 0
while n < len(names):
    x = 0
    y = 0
    z = 0
    with open(files[n]) as infile:
        for line in infile:
            splitline = line.split(",")
            for entry in splitline:
                grids[n][x][y][z] = float(entry)
                z = z + 1
                if z == 60:
                    z = 0
                    y = y + 1
                    if y == 60:
                        y = 0
                        x = x + 1
    n = n + 1

n = -1
difference_grids = []
difference_grid_names = []
a = 0
while a < len(names):
    b = 0
    while b < len(names):
        if a < b:
            difference_grids.append([])
            if lock_z_axis:
                difference_grid_names.append("Difference in stability")
            else:
                difference_grid_names.append(f"Difference in stability: {names[b]} - {names[a]}")
            n = n + 1
            x = 0
            while x < pixel_count:
                difference_grids[n].append([])
                x = x + 1
            x = 0
            while x < pixel_count:
                y = 0
                while y < pixel_count:
                    difference_grids[n][x].append([])
                    y = y + 1
                x = x + 1
            x = 0
            while x < pixel_count:
                y = 0
                while y < pixel_count:
                    z = 0
                    while z < pixel_count:
                        difference_grids[n][x][y].append(grids[b][x][y][z] - grids[a][x][y][z])
                        z = z + 1
                    y = y + 1
                x = x + 1
        b = b + 1
    a = a + 1

if plot:
    plotaxisindex = 0
    while plotaxisindex < pixel_count:
        plotvalues = []
        index1 = 0
        while index1 < pixel_count:
            plotvalues.append([])
            index1 = index1 + 1
        index1 = 0
        while index1 < pixel_count:
            index2 = 0
            while index2 < pixel_count:
                plotvalues[index1].append(0)
                if plotaxis == 0:
                    plotvalues[index1][index2] = difference_grids[plot_grid][plotaxisindex][index1][index2]
                elif plotaxis == 1:
                    plotvalues[index1][index2] = difference_grids[plot_grid][index2][plotaxisindex][index1]
                else:
                    plotvalues[index1][index2] = difference_grids[plot_grid][index1][index2][plotaxisindex]
                index2 = index2 + 1
            index1 = index1 + 1
        plt.figure()
        if lock_z_axis:
            seaborn.heatmap(plotvalues, vmin = -0.1, vmax = 0.5)
        else:
            seaborn.heatmap(plotvalues)
        plt.title(f"{Variables[plotaxis]} = index {plotaxisindex}")
        plt.ylabel(Variables[(plotaxis + 1) % 3])
        plt.xlabel(Variables[(plotaxis + 2) % 3])
    
        plotaxisindex = plotaxisindex + 1
        print(f"Progress: processed {plotaxisindex} out of {pixel_count} slices.")

if summary_one_cube:
    means_along_x = []
    medians_along_x = []
    percentiles25 = []
    percentiles75 = []
    stdevs_along_x = []
    all_values = []
    diagonal_values = []
    x = 0
    while x < pixel_count:
        all_in_x = []
        y = 0
        while y < pixel_count:
            z = 0
            while z < pixel_count:
                all_values.append(difference_grids[summary_which_cube][x][y][z])
                all_in_x.append(difference_grids[summary_which_cube][x][y][z])
                if (x == y) and (y == z):
                    diagonal_values.append(difference_grids[summary_which_cube][x][y][z])
                z = z + 1
            y = y + 1
        mean_in_x = 0
        median_in_x = []
        values_in_x = []
        for entry in all_in_x:
            mean_in_x = mean_in_x + entry
            median_in_x.append(entry)
            values_in_x.append(entry)
        stdev_in_x = np.std(values_in_x)
        mean_in_x = mean_in_x / len(all_in_x)
        median_in_x.sort()
        median_in_x_value = statistics.median(median_in_x)
        means_along_x.append(mean_in_x)
        medians_along_x.append(median_in_x_value)
        percentiles25.append(statistics.median(median_in_x[:int(len(median_in_x)/2)]))
        percentiles75.append(statistics.median(median_in_x[int(len(median_in_x)/2):]))
        median_in_x.clear()
        stdevs_along_x.append(stdev_in_x)
        x = x + 1
    mean_value = 0
    negative_values = 0
    positive_values = 0
    for entry in all_values:
        mean_value = mean_value + entry
        if entry > 0:
            positive_values = positive_values + 1
        if entry < 0:
            negative_values = negative_values + 1
    mean_value = mean_value / len(all_values)
    print(f"Mean {difference_grid_names[summary_which_cube]} = {mean_value}")
    print(sci.wilcoxon(all_values))
    print(f"Positive values = {positive_values}")
    print(f"Negative values = {negative_values}")
    
    all_values.sort()
    if len(all_values) % 2 == 1:
        median_value = all_values[len(all_values/2-0.5)]
        print(f"Median {difference_grid_names[summary_which_cube]} = {median_value}")
    else:
        median_value = (all_values[int(len(all_values)/2)] + all_values[int(len(all_values)/2-1)])/2
        print(f"Median {difference_grid_names[summary_which_cube]} = {median_value}")
    
    plt.figure()
    plt.plot(diagonal_values)
    plt.xlabel("Relative abundance of microbial species (index)")
    plt.ylabel(f"{difference_grid_names[summary_which_cube]}")
    plt.grid(axis='y')
    #plt.ylim(-0.25,0.25)
    plt.show()
    
    Mins = []
    Maxs = []
    Sems = []
    for DOC_value in difference_grids[0]:
        DOCValues = []
        Stab_min = 10
        Stab_max = -10
        for Aut_value in DOC_value:
            for Het_value in Aut_value:
                DOCValues.append(Het_value)
                if Het_value < Stab_min:
                    Stab_min = Het_value
                if Het_value > Stab_max:
                    Stab_max = Het_value
        Sems.append(sci.sem(DOCValues))
        Mins.append(Stab_min)
        Maxs.append(Stab_max)
    
    index = 0
    Sems_lower = []
    Sems_upper = []
    for entry in means_along_x:
        Sems_lower.append(means_along_x[index] - 2*Sems[index])
        Sems_upper.append(means_along_x[index] + 2*Sems[index])
        index = index + 1
    
    plt.figure()
    plt.semilogx(DOC_intermediate_biomasses, means_along_x)
    #plt.fill_between(DOC_intermediate_biomasses, Sems_lower, Sems_upper, alpha=0.3)
    plt.title("b)",loc="left")
    plt.xlabel("Detritus mass (grams per square metre)")
    plt.ylabel("Mean destabilisation effect")
    plt.grid(axis='y')
    plt.show()
    
    plt.figure()
    plt.hist(all_values, range=(-0.4,1.2), bins=32)
    plt.title("a)",loc="left")
    plt.xlabel("Difference in instability (s)")
    plt.ylabel("Frequency")
    plt.grid(True, linestyle="--")
    plt.show()
    
    negative_values = 0
    positive_values = 0
    for entry in all_values:
        if entry < 0:
            negative_values = negative_values + 1
        if entry > 0:
            positive_values = positive_values + 1
    print(f"Negative values: {negative_values}; Positive values: {positive_values}; Total values: {len(all_values)}")

if output_detritus_curve:
    with open(outputfile,"a+") as outfile:
        index = 0
        outfile.write("Mean stability,Standard deviation in stability,DOC,Median stability,25% of stability,75% of stability\n")
        while index < len(DOC_intermediate_biomasses):
            outfile.write(str(means_along_x[index]) + "," + str(stdevs_along_x[index]) + "," + str(DOC_intermediate_biomasses[index]) + "," + str(medians_along_x[index]) + "," + str(percentiles25[index]) + "," + str(percentiles75[index]) + "\n")
            index = index + 1