# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 10:38:54 2024

@author: enf21

This code is associated with, "Sloppy feeding by predators destabilises early animal food webs in the Ediacaran".
"""

import math
import seaborn
import matplotlib.pyplot as plt

write = True
outputfile = "./processed.txt"

plot = True
plotaxis = 0

plotLoopIDs = False # Do we plot loop ID variation in the cube?

loopsummary = True
PositiveLoops = {}
NegativeLoops = {}

inputfile = ".output.txt"

Variables = ["DOC biomass", "Autotrophic plankton biomass", "Heterotrophic plankton biomass"]

Variable1 = Variables[0]
Variable2 = Variables[1]
Variable3 = Variables[2]

cells = 60

entries = []
counts = []
V1min = -1
V2min = -1
V3min = -1
V1max = -1
V2max = -1
V3max = -1

var1Values = set()
var2Values = set()
var3Values = set()

if plotLoopIDs:
    PLoopIDs = []
    NLoopIDs = []

dconstindexA = 0
while dconstindexA < cells:
    entries.append([])
    counts.append([])
    if plotLoopIDs:
        PLoopIDs.append([])
        NLoopIDs.append([])
    dconstindexB = 0
    while dconstindexB < cells:
        entries[dconstindexA].append([])
        counts[dconstindexA].append([])
        if plotLoopIDs:
            PLoopIDs[dconstindexA].append([])
            NLoopIDs[dconstindexA].append([])
        dconstindexC = 0
        while dconstindexC < cells:
            entries[dconstindexA][dconstindexB].append(0)
            counts[dconstindexA][dconstindexB].append(0)
            if plotLoopIDs:
                PLoopIDs[dconstindexA][dconstindexB].append(-1)
                NLoopIDs[dconstindexA][dconstindexB].append(-1)
            dconstindexC = dconstindexC + 1
        dconstindexB = dconstindexB + 1
    dconstindexA = dconstindexA + 1

First = True
with open(inputfile) as infile:
    for line in infile:
        index = 0
        splitline = line.split(",")
        if First:
            for entry in splitline:
                if entry.strip() == Variable1:
                    Var1Column = index
                elif entry.strip() == Variable2:
                    Var2Column = index
                elif entry.strip() == Variable3:
                    Var3Column = index
                elif entry.strip() == "Stability":
                    StabilityColumn = index
                elif entry.strip() == "Heaviest positive 3 loop":
                    PositiveLoopColumn = index
                elif entry.strip() == "Heaviest negative 3 loop":
                    NegativeLoopColumn = index
                index = index + 1
            First = False
        else:
            for entry in splitline:
                if index == Var1Column:
                    Var1Value = float(entry)
                    var1Values.add(Var1Value)
                elif index == Var2Column:
                    Var2Value = float(entry)
                    var2Values.add(Var2Value)
                elif index == Var3Column:
                    Var3Value = float(entry)
                    var3Values.add(Var3Value)
                elif index == StabilityColumn:
                    StabilityValue = float(entry)
                elif index == PositiveLoopColumn:
                    entrysplit = entry.split("-")
                    newentry = []
                    for split in entrysplit:
                        newentry.append(int(split))
                    if newentry[0] == max(newentry):
                        entry = entry
                    elif newentry[1] == max(newentry):
                        entry = str(newentry[1]) + "-" + str(newentry[2]) + "-" + str(newentry[0])
                    else:
                        entry = str(newentry[2]) + "-" + str(newentry[0]) + "-" + str(newentry[1])
                    if entry in PositiveLoops:
                        PositiveLoops.update({entry:(PositiveLoops.get(entry) + 1)})
                    else:
                        PositiveLoops.update({entry:1})
                elif index == NegativeLoopColumn:
                    entrysplit = entry.split("-")
                    newentry = []
                    for split in entrysplit:
                        newentry.append(int(split))
                    if newentry[0] == max(newentry):
                        entry = entry
                    elif newentry[1] == max(newentry):
                        entry = str(newentry[1]) + "-" + str(newentry[2]) + "-" + str(newentry[0])
                    else:
                        entry = str(newentry[2]) + "-" + str(newentry[0]) + "-" + str(newentry[1])
                    if entry in NegativeLoops:
                        NegativeLoops.update({entry:(NegativeLoops.get(entry) + 1)})
                    else:
                        NegativeLoops.update({entry:1})
                index = index + 1
            if (Var1Value < V1min) or (V1min == -1):
                V1min = Var1Value
            if (Var1Value > V1max) or (V1max == -1):
                V1max = Var1Value
            if (Var2Value < V2min) or (V2min == -1):
                V2min = Var2Value
            if (Var2Value > V2max) or (V2max == -1):
                V2max = Var2Value
            if (Var3Value < V3min) or (V3min == -1):
                V3min = Var3Value
            if (Var3Value > V3max) or (V3max == -1):
                V3max = Var3Value
    
# Hardcoded values to allow for systematic sweep.
V1max = 24000
V1min = 2.4
V2max = 2100 #800
V2min = 0.21 #0.08
V3max = 1100 #400
V3min = 0.11 #0.04

V1Range = math.log(V1max) - math.log(V1min)
V2Range = math.log(V2max) - math.log(V2min)
V3Range = math.log(V3max) - math.log(V3min)
V1Step = V1Range / cells
V2Step = V2Range / cells
V3Step = V3Range / cells

boundaries1 = []
boundaries2 = []
boundaries3 = []

index = 0
while index < cells:
    boundaries1.append(math.exp(math.log(V1min) + index * V1Step))
    boundaries2.append(math.exp(math.log(V2min) + index * V2Step))
    boundaries3.append(math.exp(math.log(V3min) + index * V3Step))
    index = index + 1
boundaries1.append(V1max)
boundaries2.append(V2max)
boundaries3.append(V3max)

if plotLoopIDs:
    nPositiveLoopIDs = len(list(PositiveLoops))
    nNegativeLoopIDs = len(list(NegativeLoops))

with open(inputfile) as infile:
    First = True
    for line in infile:
        if First:
            First = False
        else:
            splitline = line.split(",")
            index = 0
            for entry in splitline:
                if index == Var1Column:
                    Var1Value = float(entry)
                elif index == Var2Column:
                    Var2Value = float(entry)
                elif index == Var3Column:
                    Var3Value = float(entry)
                elif index == StabilityColumn:
                    StabilityValue = abs(float(entry))
                elif index == PositiveLoopColumn:
                    entrysplit = entry.split("-")
                    newentry = []
                    for split in entrysplit:
                        newentry.append(int(split))
                    if newentry[0] == max(newentry):
                        entry = entry
                    elif newentry[1] == max(newentry):
                        entry = str(newentry[1]) + "-" + str(newentry[2]) + "-" + str(newentry[0])
                    else:
                        entry = str(newentry[2]) + "-" + str(newentry[0]) + "-" + str(newentry[1])
                    PositiveLoop = entry
                elif index == NegativeLoopColumn:
                    entrysplit = entry.split("-")
                    newentry = []
                    for split in entrysplit:
                        newentry.append(int(split))
                    if newentry[0] == max(newentry):
                        entry = entry
                    elif newentry[1] == max(newentry):
                        entry = str(newentry[1]) + "-" + str(newentry[2]) + "-" + str(newentry[0])
                    else:
                        entry = str(newentry[2]) + "-" + str(newentry[0]) + "-" + str(newentry[1])
                    NegativeLoop = entry
                index = index + 1
            index1 = 0
            index2 = 0
            index3 = 0
            while Var1Value > boundaries1[index1]:
                index1 = index1 + 1
            while Var2Value > boundaries2[index2]:
                index2 = index2 + 1
            while Var3Value > boundaries3[index3]:
                index3 = index3 + 1
            entries[index1-1][index2-1][index3-1] = entries[index1-1][index2-1][index3-1] + StabilityValue
            counts[index1-1][index2-1][index3-1] = counts[index1-1][index2-1][index3-1] + 1
            if plotLoopIDs:
                PLoopIDs[index1-1][index2-1][index3-1] = list(PositiveLoops).index(PositiveLoop)
                NLoopIDs[index1-1][index2-1][index3-1] = list(NegativeLoops).index(NegativeLoop)

index1 = 0
while index1 < cells:
    index2 = 0
    while index2 < cells:
        index3 = 0
        while index3 < cells:
            if counts[index1][index2][index3] > 0:
                entries[index1][index2][index3] = entries[index1][index2][index3] / counts[index1][index2][index3]
            index3 = index3 + 1
        index2 = index2 + 1
    index1 = index1 + 1

if write:
    First = True
    with open(outputfile,"w+") as outfile:
        for square in entries:
            for line in square:
                for value in line:
                    if First:
                        outfile.write(str(value))
                        First = False
                    else:
                        outfile.write(",")
                        outfile.write(str(value))

if loopsummary:
    print("Positive loop,Frequency")
    for value in PositiveLoops:
        print(f"{value},{PositiveLoops.get(value)}")
    print("")
    print("Negative loop,Frequency")
    for value in NegativeLoops:
        print(f"{value},{NegativeLoops.get(value)}")

Variables = ["Detritus biomass", "Autotrophic bacterioplankton biomass", "Heterotrophic bacterioplankton biomass"]

if plot:
    plotaxisindex = 0
    while plotaxisindex < cells:
        plotvalues = []
        index1 = 0
        while index1 < cells:
            plotvalues.append([])
            index1 = index1 + 1
        index1 = 0
        while index1 < cells:
            index2 = 0
            while index2 < cells:
                plotvalues[index1].append(0)
                if plotaxis == 0:
                    plotvalues[index1][index2] = entries[plotaxisindex][index1][index2]
                elif plotaxis == 1:
                    plotvalues[index1][index2] = entries[index2][plotaxisindex][index1]
                else:
                    plotvalues[index1][index2] = entries[index1][index2][plotaxisindex]
                index2 = index2 + 1
            index1 = index1 + 1
        plt.figure()
        #seaborn.heatmap(plotvalues, norm=LogNorm(vmin = 0.00001, vmax = 1))
        seaborn.heatmap(plotvalues, vmin = 0.00001, vmax = 1)
        plt.ylabel(Variables[(plotaxis + 1) % 3])
        plt.xlabel(Variables[(plotaxis + 2) % 3])
        if plotaxis == 0:
            plt.title(f"{Variables[plotaxis]} = {round((boundaries1[plotaxisindex] + boundaries1[plotaxisindex + 1]) / 2,2)}")
        elif plotaxis == 1:
            plt.title(f"{Variables[plotaxis]} = {(boundaries2[plotaxisindex] + boundaries2[plotaxisindex + 1]) / 2}")
        else:
            plt.title(f"{Variables[plotaxis]} = {(boundaries3[plotaxisindex] + boundaries3[plotaxisindex + 1]) / 2}")
        plt.xticks([0,10,20,30,40,50,60], labels = ["0.12","0.47","2.20","10.2","47.4","220","1021"])
        plt.yticks([0,10,20,30,40,50,60], labels = ["0.23","0.91","4.20","19.5","90.5","420","1951"])
    
        plotaxisindex = plotaxisindex + 1
        print(f"Progress: processed {plotaxisindex} out of {cells} slices.")

if plotLoopIDs:
    plotaxisindex = 0
    while plotaxisindex < cells:
        plotvalues = []
        index1 = 0
        while index1 < cells:
            plotvalues.append([])
            index1 = index1 + 1
        index1 = 0
        while index1 < cells:
            index2 = 0
            while index2 < cells:
                plotvalues[index1].append(0)
                if plotaxis == 0:
                    plotvalues[index1][index2] = PLoopIDs[plotaxisindex][index1][index2]
                elif plotaxis == 1:
                    plotvalues[index1][index2] = PLoopIDs[index2][plotaxisindex][index1]
                else:
                    plotvalues[index1][index2] = PLoopIDs[index1][index2][plotaxisindex]
                index2 = index2 + 1
            index1 = index1 + 1
        plt.figure()
        #seaborn.heatmap(plotvalues, norm=LogNorm(vmin = 0.00001, vmax = 1))
        seaborn.heatmap(plotvalues, vmin = -1, vmax = len(list(PositiveLoops)))
        if plotaxis == 0:
            plt.title(f"Sliced plot, {Variables[plotaxis]} = {(boundaries1[plotaxisindex] + boundaries1[plotaxisindex + 1]) / 2}")
        elif plotaxis == 1:
            plt.title(f"Sliced plot, {Variables[plotaxis]} = {(boundaries2[plotaxisindex] + boundaries2[plotaxisindex + 1]) / 2}")
        else:
            plt.title(f"Sliced plot, {Variables[plotaxis]} = {(boundaries3[plotaxisindex] + boundaries3[plotaxisindex + 1]) / 2}")
        plt.ylabel(Variables[(plotaxis + 1) % 3])
        plt.xlabel(Variables[(plotaxis + 2) % 3])
    
        plotaxisindex = plotaxisindex + 1
        print(f"Progress: processed {plotaxisindex} out of {cells} slices.")
        
    plotaxisindex = 0
    while plotaxisindex < cells:
        plotvalues = []
        index1 = 0
        while index1 < cells:
            plotvalues.append([])
            index1 = index1 + 1
        index1 = 0
        while index1 < cells:
            index2 = 0
            while index2 < cells:
                plotvalues[index1].append(0)
                if plotaxis == 0:
                    plotvalues[index1][index2] = NLoopIDs[plotaxisindex][index1][index2]
                elif plotaxis == 1:
                    plotvalues[index1][index2] = NLoopIDs[index2][plotaxisindex][index1]
                else:
                    plotvalues[index1][index2] = NLoopIDs[index1][index2][plotaxisindex]
                index2 = index2 + 1
            index1 = index1 + 1
        plt.figure()
        #seaborn.heatmap(plotvalues, norm=LogNorm(vmin = 0.00001, vmax = 1))
        seaborn.heatmap(plotvalues, vmin = -1, vmax = len(list(NegativeLoops)))
        if plotaxis == 0:
            plt.title(f"Sliced plot, {Variables[plotaxis]} = {(boundaries1[plotaxisindex] + boundaries1[plotaxisindex + 1]) / 2}")
        elif plotaxis == 1:
            plt.title(f"Sliced plot, {Variables[plotaxis]} = {(boundaries2[plotaxisindex] + boundaries2[plotaxisindex + 1]) / 2}")
        else:
            plt.title(f"Sliced plot, {Variables[plotaxis]} = {(boundaries3[plotaxisindex] + boundaries3[plotaxisindex + 1]) / 2}")
        plt.ylabel(Variables[(plotaxis + 1) % 3])
        plt.xlabel(Variables[(plotaxis + 2) % 3])
    
        plotaxisindex = plotaxisindex + 1
        print(f"Progress: processed {plotaxisindex} out of {cells} slices.")

var1ValuesList = list(var1Values)
var1ValuesList.sort()
var2ValuesList = list(var2Values)
var2ValuesList.sort()
var3ValuesList = list(var3Values)
var3ValuesList.sort()

if plot:
    StabilityAlongDOC = []
    CountsAlongDOC = []
    StabilityAlongAut = []
    CountsAlongAut = []
    StabilityAlongHet = []
    CountsAlongHet = []
    index0 = 0
    index1 = 0
    index2 = 0
    while index0 < cells:
        StabilityAlongDOC.append(0)
        CountsAlongDOC.append(0)
        StabilityAlongAut.append(0)
        CountsAlongAut.append(0)
        StabilityAlongHet.append(0)
        CountsAlongHet.append(0)
        index0 += 1
    index0 = 0
    for DOC in entries:
        index1 = 0
        for Aut in DOC:
            index2 = 0
            for Het in Aut:
                StabilityAlongDOC[index0] = StabilityAlongDOC[index0] + Het
                CountsAlongDOC[index0] = CountsAlongDOC[index0] + 1
                StabilityAlongAut[index1] = StabilityAlongAut[index1] + Het
                CountsAlongAut[index1] = CountsAlongAut[index1] + 1
                StabilityAlongHet[index2] = StabilityAlongHet[index2] + Het
                CountsAlongHet[index2] = CountsAlongHet[index2] + 1
                index2 += 1
            index1 += 1
        index0 += 1
    index0 = 0
    while index0 < cells:
        StabilityAlongDOC[index0] = StabilityAlongDOC[index0] / CountsAlongDOC[index0]
        StabilityAlongAut[index0] = StabilityAlongAut[index0] / CountsAlongAut[index0]
        StabilityAlongHet[index0] = StabilityAlongHet[index0] / CountsAlongHet[index0]
        index0 += 1
    plt.figure()
    plt.semilogx(var1ValuesList, StabilityAlongDOC)
    plt.title("b)",loc="left")
    plt.xlabel("Detritus mass (grams per square metre)")
    plt.ylabel("Mean stability (s)")
    plt.grid(axis='y')
    plt.show()
    plt.figure()
    plt.semilogx(var2ValuesList, StabilityAlongAut)
    plt.title("b)",loc="left")
    plt.xlabel("Aut mass (grams per square metre)")
    plt.ylabel("Mean stability (s)")
    plt.grid(axis='y')
    plt.show()
    plt.figure()
    plt.semilogx(var3ValuesList, StabilityAlongHet)
    plt.title("b)",loc="left")
    plt.xlabel("Het mass (grams per square metre)")
    plt.ylabel("Mean stability (s)")
    plt.grid(axis='y')
    plt.show()
    
    all_values = []
    for entryA in entries:
        for entryB in entryA:
            for entryC in entryB:
                all_values.append(entryC)
    
    plt.figure()
    plt.hist(all_values, range=(0,1), bins=40)
    plt.title("a)",loc="left")
    plt.xlabel("Stability (s)")
    plt.ylabel("Frequency")
    plt.grid(True, linestyle="--")
    plt.show()
    