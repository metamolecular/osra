#!/usr/bin/python
import sys
import operator
from operator import itemgetter
from os import listdir
from os.path import isfile, join
import random
import math

from openbabel import *

#verify_model = [-0.02285364975052934, 0.20014736496836238, 0.24359177367611942, 0.09879557707190285, 0.23658236605437208, 0.07640365396444981, -0.016708707126055426, 0.291372680023203, 0.19743466243758914, -0.04143248625528295, 0.13283691497820155, -0.09404435905499846, -0.34011551678018254, -0.036998642457270414, 0.3366565862477758, 0.2528949626509886, 0.3523503866589532, 0.3013887466571869, 0.2457272087651062, -0.08224552150295372, 0.0386321456632419, 0.2269247796030229, 0.191691888047917, 0.029364205782588967, -0.09207024117341821, 0.024551588143053422]

verify_model = [-0.11469143725730054, 0.15723547931889853, 0.19765680222250673, 0.249101590474403, 0.1897669087341134, 0.19588348907301223, 0.3354622208036507, 0.16779269801176255, -0.21232000222198893, 0.016958281784354032, -0.08672059360133752, -0.05105752296619957, -0.349912750824004, 0.18836317536530647, 0.22316782354758827, 0.27741998968081166, 0.25710999274481955, 0.27968899280120096, 0.12695166847876285, -0.10020778884718293, 0.05150631410596443, 0.22283571763712148, 0.23130179826714167, 0.1049054095759948, 0.05333970810460394, -0.12491056666737535]               

#verify_model = [-0.1545855719726278, -0.16679291864636722, 0.2779073764931343, -0.183848833684335, 0.010790075194024773, 0.29094165316568404, 0.055324497605819506, -0.2104820104189514, -0.1781856691338483, 0.12170164214195042, -0.03319968208305941, -0.17050232311223057, -0.3855170942775288, -0.07088710430614285, 0.24005317771967355, 0.2759926472483148, 0.25348276233777095, 0.23427258354038655, -0.1175747967837222, -0.18681840394577787, 0.06103578120099978, 0.24422743725717977, 0.25207495568639754, -0.09625789745569688, -0.01025153552468599, 0.19182292957981223]

#verify_model = [-0.18431997080588122, -0.13632503439995766, 0.2891372503939111, 0.169268288671698, 0.07457361791041998, 0.1896239880096002, 0.17921064798323905, 0.24807741146917148, -0.1415236210886208, 0.017500171104361622, -0.1444618582502517, -0.019936025471384962, -0.39685240156986173, -0.22908080638789957, 0.27059782240339336, 0.17386007711539425, 0.21106985232185135, 0.2865651377997317, 0.09715634915474097, 0.008962730235627716, 0.030950650868271857, 0.2256707621011711, 0.19237430308515277, -0.22938527889531524, 0.15229124226660562, 0.2099925925427031]

def normalize(x):
    n = sum(map( operator.mul, x, x))
    n = math.sqrt(n)
    y = [a/n for a in x]
    return y

def trial_confidence(x,c):
    return sum(map( operator.mul, x, c))

def model_recall(N,res_iter_all,probabilities,target,inchi_list,total):
    recall_model = 0
    k = 0
    for i in range(N):
        total_probabilities = [0.,0.,0.,0.,0.]
        n_probabilities =     [0.,0.,0.,0.,0.]
        recall_inchi = set()
        for r in range(0,5):
            for j in range(len(res_iter_all[i])):
                if (res_iter_all[i][j] == r):
                    total_probabilities[r] += probabilities[k+j]
                    n_probabilities[r] += 1
        maxp = 0
        maxr = 0
        for r in range(0,5):
            if (n_probabilities[r]>0 and maxp < total_probabilities[r]/n_probabilities[r]):
                maxp = total_probabilities[r]/n_probabilities[r]
                maxr = r
        first = True
        for r in range(0,5):
            if (n_probabilities[r]>0 and maxp == total_probabilities[r]/n_probabilities[r] and (r == 2 or r == 3) and first):
                maxr = r
                first = False

        for j in range(len(res_iter_all[i])):
            if (res_iter_all[i][j] == maxr and target[k+j] == 1):
                recall_inchi.add(inchi_list[k+j])
        
        k += len(res_iter_all[i])
        recall_model += len(recall_inchi);

    return 1.*recall_model/total

def mutation(n):
    d = []
    for i in range(n):
        d.append(2.*random.random()-1.)
    return normalize(d)

def crossover(population):
    n = len(population)
    i = random.randint(0,n-1)
    j = random.randint(0,n-1)
    v = population[i][1]
    u = population[j][1]
    m = random.randint(0,len(v)-1)
    r = v[:m]
    r.extend(u[m:])
    return normalize(r)

obconversion1 = OBConversion()
obconversion1.SetInFormat("sdf")
obconversion1.SetOutFormat("inchi")
obmol1 = OBMol()
obconversion2 = OBConversion()
obconversion2.SetInFormat("sdf")
obconversion2.SetOutFormat("inchi")
obmol2 = OBMol()

result = OBPlugin.ListAsString("fingerprints")
assert "FP2" in result, result
fingerprinter = OBFingerprint.FindFingerprint("FP2")
v1 = vectorUnsignedInt()
v2 = vectorUnsignedInt()
obErrorLog.StopLogging()

path1 = sys.argv[1]
path2 = sys.argv[2]
files = [ f for f in listdir(path1) if isfile(join(path1,f)) ]

total = 0;
recall = 0;
target = []
train = []
confidence = []
resolutions = []
single = []
res_iter_all = []
probabilities = []
inchi_list = []
for f in files:
    file1 = join(path1,f);
    inchi_set1 = set()
    notatend1 = obconversion1.ReadFile(obmol1,file1)
    while notatend1:
        obmol1.AddHydrogens()
        inchi1 =  obconversion1.WriteString(obmol1)
        if inchi1:
            inchi_set1.add(inchi1)
            total += 1;
        obmol1 = OBMol()
        notatend1 = obconversion1.Read(obmol1)

    file2 = join(path2,f);
    inchi_set2 = set()
    if isfile(file2):
        notatend2 = obconversion2.ReadFile(obmol2,file2)
        data_file = []
        resolution = []
        res_iter = []
        while notatend2:
            obmol2.AddHydrogens()
            line = obmol2.GetData("Confidence_parameters").GetValue()
            data = [int(d) for d in line.split(",")]
            data_file.append(data)
            inchi2 =  obconversion2.WriteString(obmol2)
            result = 0 
            if inchi2:
                inchi_set2.add(inchi2)
            if inchi2 in inchi_set1:
                result = 1
            target.append(result)
            inchi_list.append(inchi2);
            resolution.append(int(obmol2.GetData("Resolution").GetValue()))
            res_iter.append(int(obmol2.GetData("Resolution_iteration").GetValue()))
            obmol2 = OBMol()
            notatend2 = obconversion2.Read(obmol2)
        resolutions.append(resolution)
        res_iter_all.append(res_iter)
        for d in data_file:
            train.append(d)
        recall += len(inchi_set1.intersection(inchi_set2))

N = len(files);
ideal_recall = 1.*recall/total
print "Ideal: ",ideal_recall
population = []

if len(sys.argv)>3 and sys.argv[3] == "-verify":
    c = verify_model
    probabilities = []
    for t in train:
        probabilities.append(trial_confidence(t,c))
    r = model_recall(N,res_iter_all,probabilities,target,inchi_list,total)
    print "Model: ",r
    exit(0)

for i in range(100):
    c = mutation(len(train[0]))
    probabilities = []
    for t in train:
        probabilities.append(trial_confidence(t,c))
    r = model_recall(N,res_iter_all,probabilities,target,inchi_list,total)
    population.append([r,c])
    
population.sort(key=itemgetter(0),reverse=True)
round = 1
print round,population[0][0]

while population[0][0] < ideal_recall - 0.001 and round < 1000:
    # keep top 10
    new_population = population[0:10]
    # mutation 10
    for i in range(10):
        c = mutation(len(train[0]))
        probabilities = []
        for t in train:
            probabilities.append(trial_confidence(t,c))
        r = model_recall(N,res_iter_all,probabilities,target,inchi_list,total)
        new_population.append([r,c])
    # crossover 80
    for i in range(80):
        c = crossover(population)
        probabilities = []
        for t in train:
            probabilities.append(trial_confidence(t,c))
        r = model_recall(N,res_iter_all,probabilities,target,inchi_list,total)
        new_population.append([r,c])

    new_population.sort(key=itemgetter(0),reverse=True)
    population = new_population
    round += 1
    print round,population[0][0]

print population[0][1]











