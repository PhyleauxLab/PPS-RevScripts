from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import SeqIO
import string
import shutil
import random
import glob
import os




### generate empty dataset for use in mrbayes
def generateEmptyData(numBases, numTaxa):

    alphabet = string.ascii_uppercase
    sequences = {}
    if numTaxa < len(alphabet):
        for i in range(numTaxa):
            sequences[alphabet[i]+"1"] = "N" * numBases
    
    
    ## build sequence set in NEXUS format
    records = []
    for key, val in sequences.items():
        record = SeqRecord(Seq(val, IUPAC.ambiguous_dna), id=key, name=key, description="")
        records.append(record)

    with open("pps_empty_data.nex", "w") as output_handle:
        SeqIO.write(records, output_handle, "nexus")

### run empty dataset in MrBayes
def runEmptyData():
    
    ## setup a basic folder structure
    os.mkdir("empty_data")
    os.mkdir("empty_data/gtr")
    os.mkdir("empty_data/jc")
    
    ## copy empty_data and add desired Bayesblock to the file for GTR
    gtrBayesBlock = "\nbegin mrbayes;\n\tlset nst=6;\n\tprset revratepr = symdir(1.0) statefreqpr = Dirichlet(1.0);\n\tmcmcp nchains=2 ngen=100000 samplefreq=100 printfreq=1000;\n\tmcmc;\nend;"
    shutil.copyfile("pps_empty_data.nex", "empty_data/gtr/pps_empty_data_gtr.nex")
    handle =  open("empty_data/gtr/pps_empty_data_gtr.nex", "a")
    handle.write(gtrBayesBlock)
    handle.close()         

    ## run MrBayes analysis on empty data for GTR
    os.system("mb empty_data/gtr/pps_empty_data_gtr.nex")

    ## copy empty_data and add desired Bayesblock to the file for JC
    jcBayesBlock = "\nbegin mrbayes;\n\tlset nst=1;\n\tprset statefreqpr=fixed(equal);\n\tmcmcp nchains=2 ngen=100000 samplefreq=100 printfreq=1000;\n\tmcmc;\nend;"
    shutil.copyfile("pps_empty_data.nex", "empty_data/jc/pps_empty_data_jc.nex")
    handle =  open("empty_data/jc/pps_empty_data_jc.nex", "a")
    handle.write(jcBayesBlock)
    handle.close()         

    ## run MrBayes analysis on empty data for JC
    os.system("mb empty_data/jc/pps_empty_data_jc.nex")


## Randomly Sample from results to generate 100 trees and param sets
def sampleEmptyPosteriors(burnin, numDatasets, numTaxa):
    
    ## directory structure to keep things nice and tidy
    os.mkdir("empty-data-posterior-samples")
    os.mkdir("empty-data-posterior-samples/params")
    os.mkdir("empty-data-posterior-samples/trees")
    os.mkdir("empty-data-posterior-samples/trees/gtr")
    os.mkdir("empty-data-posterior-samples/trees/jc")
    
    ## read in the gtr trees and parameters from empty data
    gtrTrees = []
    for i in range(2):
        name = "empty_data/gtr/pps_empty_data_gtr.nex.run" + str(i+1) + ".t"
        lines = open(name, "r").readlines()
        treeLines = lines[5 + numTaxa:]
        for line in treeLines[burnin:]:
            if "tree gen." in line:
                gtrTrees.append(line.strip() + "|" + str(i+1))
    
    ## read in the gtr trees and parameters from empty data
    gtrParams = []
    for i in range(2):
        name = "empty_data/gtr/pps_empty_data_gtr.nex.run" + str(i+1) + ".p"
        lines = open(name, "r").readlines()
        paramLines = lines[2:]
        for line in paramLines[burnin:]:
            if line:
                gtrParams.append(line.strip() + "\t" + str(i+1))
    
    ## randomly sample GTR trees and params and save values to new file

    ## get random values for list indices to sample
    samples = random.sample(range(0,len(gtrTrees)), numDatasets)
    gtrParamHandle = open("empty-data-posterior-samples/params/gtr_empty_data_random_param_samples.tsv", "w")
    gtrParamHandle.write("Gen\tLnL\tLnPr\tTL\tr(A<->C)\tr(A<->G)\tr(A<->T)\tr(C<->G)\tr(C<->T)\tr(G<->T)\tpi(A)\tpi(C)\tpi(G)\tpi(T)\tRun\n")
    gens = []
    for sample in samples:
        ## write the tree to a file for seq-gen
        tree = gtrTrees[sample].strip().split("= [&U]")[-1].split("|")[0]
        gen = gtrTrees[sample].strip().split("= [&U]")[0].split(".")[-1].replace(" ", "")
        run = gtrTrees[sample].strip().split("= [&U]")[-1].split("|")[-1]
        treeHandle = open("empty-data-posterior-samples/trees/gtr/gtr_random_tree_gen_"+str(gen)+"_run_"+str(run)+".tre", "w")
        treeHandle.write(tree)
        treeHandle.close()

        ## write the params to a file for seq-gen
        param = gtrParams[sample]
        gtrParamHandle.write(str(param) + "\n")

    gtrParamHandle.close()
    ## rinse and repeat for JC trees and params
    ## randomly sample JC trees and params and save values to new file

    ## read in the jc trees and parameters from empty data
    jcTrees = []
    for i in range(2):
        name = "empty_data/jc/pps_empty_data_jc.nex.run" + str(i+1) + ".t"
        lines = open(name, "r").readlines()
        treeLines = lines[5 + numTaxa:]
        for line in treeLines[burnin:]:
            if "tree gen." in line:
                jcTrees.append(line.strip() + "|" + str(i+1))
    
    ## read in the jc trees and parameters from empty data
    jcParams = []
    for i in range(2):
        name = "empty_data/jc/pps_empty_data_jc.nex.run" + str(i+1) + ".p"
        lines = open(name, "r").readlines()
        paramLines = lines[2:]
        for line in paramLines[burnin:]:
            if line:
                jcParams.append(line.strip() + "\t" + str(i+1))
    
    ## randomly sample jc trees and params and save values to new file

    ## get random values for list indices to sample
    jcParamHandle = open("empty-data-posterior-samples/params/jc_empty_data_random_param_samples.tsv", "w")
    jcParamHandle.write("Gen\tLnL\tLnPr\tTL\tRun\n")
    for sample in samples:
        ## write the tree to a file for seq-gen
        tree = jcTrees[sample].strip().split("= [&U]")[-1]
        gen = jcTrees[sample].strip().split("= [&U]")[0].split(".")[-1].replace(" ", "")
        run = gtrTrees[sample].strip().split("= [&U]")[-1].split("|")[-1]
        treeHandle = open("empty-data-posterior-samples/trees/jc/jc_random_tree_gen_"+str(gen)+"_run_"+str(run)+".tre", "w")
        treeHandle.write(tree)
        treeHandle.close()

        ## write the params to a file for seq-gen
        param = jcParams[sample]
        jcParamHandle.write(str(param) + "\n")

    jcParamHandle.close()


## Use random trees and param sets to generate sequences
def simulateSequences(numTaxa, numBases, numDatasets):
    
    ## create directories to store simulated sequences
    os.mkdir("simulated-sequences")
    os.mkdir("simulated-sequences/gtr")
    os.mkdir("simulated-sequences/jc")

    ## simulate sequence data using priors from empty data analysis under GTR
    ## pull tree and params from each file.
    gtrTreeFiles = glob.glob("empty-data-posterior-samples/trees/gtr/*.tre")
    gtrParams = {}
    gtrParamLines = open("empty-data-posterior-samples/params/gtr_empty_data_random_param_samples.tsv").readlines()
    for line in gtrParamLines[1:]:
        rec = line.strip().split()

        ## reformats the values to seq-gen friendly floats / strings
        rec4 = "{0:.8f}".format(float(rec[4]))
        rec5 = "{0:.8f}".format(float(rec[5]))
        rec6 = "{0:.8f}".format(float(rec[6]))
        rec7 = "{0:.8f}".format(float(rec[7]))
        rec8 = "{0:.8f}".format(float(rec[8]))
        rec9 = "{0:.8f}".format(float(rec[9]))

        rec10 = "{0:.8f}".format(float(rec[10]))
        rec11 = "{0:.8f}".format(float(rec[11]))
        rec12 = "{0:.8f}".format(float(rec[12]))
        rec13 = "{0:.8f}".format(float(rec[13]))
        gtrParams[rec[0]] = "-r"+str(rec4)+ " " +str(rec5)+ " " +str(rec6)+ " " +str(rec7)+ " " +str(rec8)+ " " +str(rec9)+" -f"+str(rec10)+ " " +str(rec11)+ " " +str(rec12)+ " " +str(rec13)
    
    ## call seqgen in appropriate way to simulate sequence set
    count = 1
    for tree in gtrTreeFiles:
        gen = tree.split("gen_")[-1].split("_")[0]
        run = tree.split("_")[-1].replace(".tre", "")
        seqName = "simulated-sequences/gtr/gtr_simulated_dataset_" + str(count) + "_gen_" + str(gen) + "_run_" + str(run) + ".nex"
        gtrCMD = "./seq-gen -mGTR -l" + str(numBases) + " " + gtrParams[gen] + " -on " + tree + " > " + seqName
        os.system(gtrCMD)
        count += 1
    
    ## simulate sequence data using priors from empty data analysis under JC
    ## pull tree and params from each file.
    jcTreeFiles = glob.glob("empty-data-posterior-samples/trees/jc/*.tre")
    ## call seqgen in appropriate way to simulate sequence set
    count = 1
    for tree in jcTreeFiles:
        gen = tree.split("gen_")[-1].split("_")[0]
        run = tree.split("_")[-1].replace(".tre", "")
        seqName = "simulated-sequences/jc/jc_simulated_dataset_" + str(count) + "_gen_" + str(gen) + "_run_" + str(run) + ".nex"
        jcCMD = "./seq-gen -mHKY -l" + str(numBases) + " -on " + tree + " > " + seqName
        os.system(jcCMD)
        count += 1
    

####################################
## set these variables before use ##
####################################
numBases = 1000
numDatasets = 100
numTaxa = 16
burnin = 250 # in num of trees per run


generateEmptyData(numBases, numTaxa)
runEmptyData()
sampleEmptyPosteriors(burnin, numDatasets, numTaxa)
simulateSequences(numTaxa, numBases, numDatasets)