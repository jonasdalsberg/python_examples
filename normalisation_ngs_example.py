#!/usr/bin/env python3

### Normalisation of metagenomics data from different cheeses
### For a project in a course on next generation sequencing (NGS) at DTU
### Code written by Jonas Dalsberg Jørgensen in January 2022

import math

infile = open("cheese.data.dat", "r")
annofile = open("cheese.annotations.dat", "r")
outfile = open("cheese.normalised.dat", "w")
log_outfile = open("cheese.log.normalised.dat", "w")

linecount = 0
data = dict()
sums = list()

# annotations
annodict = dict()
cnt = 0
for line in annofile:
    cnt += 1
    if cnt > 1:
        annoline = line.strip().split("\t")
        if len(annoline) == 21:
            annodict[annoline[0]] = annoline[1]

# data collection
for line in infile:
    sum = 0
    linecount += 1
    if linecount > 1:
        genename = line.strip().split("\t")[0]
        data[genename] = []
        for num in line.strip().split("\t")[1:]:
            data[genename].append(int(num))
            sum += int(num)
        sums.append(sum)

# across sample normalisation per gene (relative to the average value across samples)
cnt = -1
newdata = dict()
for genename in data:
    cnt += 1
    if sums[cnt] != 0:
        newdata[genename] = []
        for num in data[genename]:
            newdata[genename].append(num / (sums[cnt]/16))

# compensation for sample size normalisation
read_list = [6488718,8641776,4413162,1148752,8979654,6212268,8518820,3113788,10451042,2378334,8085990,5290682,11090370,7366116,4882496,8606638]
read_sum = 0
for read in read_list:
    read_sum += read
read_avg = read_sum / 16

for i in range(16):
    read_list[i] = read_list[i] / read_avg

for genename in newdata:
    for i in range(16):
        newdata[genename][i] = newdata[genename][i] / read_list[i]

# average across cheese types
finaldata = dict()
cnt = 0
sum = 0
for genename in newdata:
    finaldata[genename] = []
    for i in range(16):
        cnt += 1
        sum += newdata[genename][i]
        if cnt == 4:
            finaldata[genename].append(sum)
            cnt = 0
            sum = 0

# top 5 most abundant genes:
maxdict = dict()
maxval = 0
cnt = -1
for genename in data:
    cnt += 1
    if sums[cnt] > maxval:
        maxdict[genename] = sums[cnt]
        maxval = sums[cnt]
maxgenes = ['smok_4_k141_901_38', 'smok_3_k141_1543_127', 'blue_3_k141_3085_17', 'blue_1_k141_5954_6', 'blue_1_k141_4164_1']

# write to file (log)
outfile.write("Gene\tBlue_cheese\tCheddar\tSmoked_cheese\tTomme\n")
for genename in finaldata:
    if genename in maxgenes:
        outfile.write("{}\t{}\t{}\t{}\t{}\n".format(annodict[genename + "_0"], finaldata[genename][0], finaldata[genename][1], finaldata[genename][2], finaldata[genename][3]))

# zero replacement and log(2) transformation
for genename in finaldata:
    for i in range(4):
        if finaldata[genename][i] == 0:
            finaldata[genename][i] += 0.000000000000000000000001
        finaldata[genename][i] = math.log2(finaldata[genename][i])

# write to file (log)
log_outfile.write("Gene\tBlue_cheese\tCheddar\tSmoked_cheese\tTomme\n")
for genename in finaldata:
    if genename in maxgenes:
        log_outfile.write("{}\t{}\t{}\t{}\t{}\n".format(annodict[genename + "_0"], finaldata[genename][0], finaldata[genename][1], finaldata[genename][2], finaldata[genename][3]))

# close files
infile.close()
annofile.close()
outfile.close()
log_outfile.close()

#print(data)
#print(sums)
#print(newdata)
#print(len(data))
#print(len(newdata))
