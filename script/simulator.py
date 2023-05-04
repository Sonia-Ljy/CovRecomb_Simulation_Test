
import random
from progress.bar import Bar
from numpy.random import default_rng
import csv
import os
from numpy import *
import pandas as pd
from Bio import SeqIO
import numpy as np
os.chdir(os.getcwd()+"/")


def generate_sequence(length):
    """
    This method will generate a sequence, and set the Sequence object's 
    sequence to that sequence.
    """
    sequence = ''
    for i in range(length):
        import random
        letter = random.choice(['A', 'T', 'G', 'C'])
        sequence += letter
    return sequence


def lin_initial_sequence(pop_size,code_path,length = 1000,ref = True):
    # sequence_linA = generate_sequence(length)
    with open(code_path+"EPI_ISL_402125.fasta", "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            sequence_linA = str(record.seq)

    lin_seq = {}
    if ref == False:
        num = 0
        for i in range(pop_size):
            num += 1
            lin_ID = "lin"+ str(num)
            lin_seq[lin_ID] = sequence_linA
    elif ref == True:
        lin_seq["reference"] = sequence_linA
    return lin_seq


def creat_dir(turns_file):
    import os
    if not os.path.exists(turns_file):
        os.makedirs(turns_file)


def choose_new_letter(letter):
    """
    This function chooses a new letter from ATGC that is
    different from the letter passed into the function.

    INPUTS:
    -	CHAR: letter
            the letter that will not be chosen from ATGC.
    """
    possible_letters = set(['A', 'T', 'G', 'C'])
    import random
    new_letter = random.choice(list(
        possible_letters.difference(set(letter)))) #随机选取除了原碱基之外的一个碱基
    return new_letter


rng = default_rng() # Random number generator
def mutate(sequence, mutation_site, qmatrix):
    '''
        A function that mutates a selected individual sequence and checks
        if the mutation was synonymous or nonsynonymous. Individual fitness
        gets decreased for each nonsynonymous mutation.

        Parameters:
            individual (dict): A dictionary representing an individual in
                               a population.
            mutation_sites (list): Sites where mutations occur.
            positions (dict): Dictionary with information about each
                                  coding position in a genome.

        Returns:
            individual (dict): Returns individual dictionary with a
                               mutated sequence and changed fitness
                               value for each nonsynonymous mutation.
    '''
    # if multiple sites in a single codon, check one site at a time
    # mutation matrix to be made customizable
    site = int(mutation_site)
    mutated_base = sequence[site]
    probs = qmatrix[mutated_base]
    possibleBases = list(probs.keys())
    probBases = list(probs.values())
    new_base =  rng.choice(possibleBases, p = probBases)
    # new_base = choose_new_letter(mutated_base)
    return new_base


# Initialize sequences
def generate_ref_and_lin(seq_length, code_path,parallel_prop, seed_gen, seed_mut_perlin,seed_gene_homo,num_seeds):
    ref_seq = lin_initial_sequence(1, code_path,seq_length, ref = True)
    temp_pop = {}
    all_pos = [x for x in range(0, seq_length)]
    import random
    for gen in range(seed_gen):
        gen += 1
        homo_mut = random.sample(all_pos, seed_gene_homo)
        seq_homo = ""
        for site in all_pos:
            if site in homo_mut:
                seq_homo += choose_new_letter(ref_seq["reference"][site])
            else:
                seq_homo += ref_seq["reference"][site]

        mut_lin_num = int((2**gen)*parallel_prop)
        import random
        select_seq_with_homo = random.sample(range(2**gen), mut_lin_num)
        already_pos = []

        for n in range(2**gen):
            seed_name = "temp_gen"+str(gen)+"_seq"+str(n)
            left_pos = [x for x in all_pos if x not in already_pos]
            import random
            positions = random.sample(left_pos, seed_mut_perlin)
            already_pos.extend(positions)
            seed_seq = ""
            if gen >= 2:
                selected_seed_name = temp_pop["temp_gen"+str(int(gen-1))+"_seq"+str(int(n/2))]
            elif gen == 1:
                selected_seed_name = ref_seq["reference"]
                
            for position in all_pos:
                if position in positions:
                    letter = choose_new_letter(selected_seed_name[position])
                    seed_seq += letter
                elif position in homo_mut and n in select_seq_with_homo:
                    seed_seq += seq_homo[position]
                else:
                    letter = selected_seed_name[position]
                    seed_seq += letter
            temp_pop[seed_name] = seed_seq

    last_pop = []
    for i in temp_pop:
        if "gen"+str(seed_gen) in i :
            last_pop.append(i)
    import random
    seed_name = random.sample(last_pop, num_seeds)
    seed_pop = {}
    for i in range(1,num_seeds+1):
        seed_pop["lin"+str(i)] = temp_pop[seed_name[i-1]]
        
    return ref_seq, seed_pop


def output_fv_file(generations, lin_fv75, turns_file):
    with open(turns_file+"/fv_75.txt", "w") as f:
        for l in lin_fv75:
            f.write(l+","+str(generations)+",")
            for v in lin_fv75[l]:
                f.write(v+",")
            f.write("\n")


def extract_fv(DIRPATH):
    lineage_file = 'fv_75_norm_cluster75.txt'
    Lineage_v = {}
    linA_list = []
    feature_mutations = []
    with open(DIRPATH+lineage_file,'r') as f:
        for i in csv.reader(f):
            if len(i[2:])>= 1:
                linA_list.append(i[0])
                Lineage_v[i[0]] = i[2:]
                for v in i[2:]:
                    if v not in feature_mutations:
                        feature_mutations.append(v)
    return Lineage_v,linA_list,feature_mutations


def extract_ref_mut(ref_seq, seq_samples, lin_all_mutate, s):
    for g in range(len(seq_samples[s])):
        if seq_samples[s][g] != ref_seq["reference"][g]:
            mut = str(g+1)+"_"+str(seq_samples[s][g])
            lin_all_mutate.append(mut)
    return lin_all_mutate


def get_lin_name(x1_lin):
    if "bk" in x1_lin:
        lin_name= ""
        flag = 0
    elif x1_lin.count("_") == 0:
        lin_name = x1_lin
        flag = 1
    else:
        lin_name = "lin"+x1_lin.split("lin")[1]
        flag = 1

    return lin_name,flag


def generate_recom_bk(seq_length,generation,x1_lin,x2_lin,pop):
    breakup = random.sample(range(0, seq_length), 1)
    new_recom_name = "gen"+str(generation)+"_"+x1_lin + "bk"+str(breakup[0])+x2_lin+"/"
    # add seq:
    left = pop[x1_lin][:breakup[0]] 
    right = pop[x2_lin][breakup[0]:]
    recom_seq = left + right
    pop[new_recom_name] = recom_seq 
    return pop


def simulator_withhomo(turns_list,seed_mut_perlin,generations_list, rec_rate):
    # print('process_id: (%s)...' % ( os.getpid()))
    for genera in generations_list:
        import copy
        turns_list_deep = copy.deepcopy(turns_list)
        for turns in range(1,turns_list_deep+1):
            ## generate ref_seq and lineage_seq
            ref_seq, all_pop = generate_ref_and_lin(seq_length, code_path,parallel_prop, seed_gen, seed_mut_perlin,seed_gene_homo,num_seeds)
            ## start the No. test
            seq_samples = []
            import copy
            pop = copy.deepcopy(all_pop)
            # Loop for every process
            for generation in range(1, genera + 1):
                all_pos = [x for x in range(0,seq_length)]
                import random
                homo_mut = random.sample(all_pos, 1)
                seq_homo = ""
                for site in all_pos:
                    if site in homo_mut:
                        seq_homo += mutate(ref_seq["reference"],site,qmatrix)
                    else:
                        seq_homo += ref_seq["reference"][site]

                mut_lin_num = int(len(pop)*parallel_prop) 
                select_seq_with_homo = random.sample(range(len(pop)),mut_lin_num) 
                ## Mutation
                pool = copy.deepcopy(pop)
                count = 0
                for lin in pool:
                    count+=1
                    if "bk" not in lin:
                        new_lin_name = "gen"+str(generation)+"_"+lin
                        seed_sequence = pop[lin]
                        mutations_seq = "" # The mutated sequence based on the seed (or former) sequence
                        # Calculated the number of mutate sites and positions based on the length of genomes.
                        num_mut = round(sum(rng.poisson(mut_rate,seq_length))) # rate is per site

                        if num_mut > 0:
                            import random
                            positions = random.sample(range(0, seq_length), num_mut)

                            if (count-1) not in select_seq_with_homo:
                                for position in range(0,len(seed_sequence)):
                                    if position in positions:
                                        letter = mutate(seed_sequence, position,qmatrix)
                                        mutations_seq += letter
                                    else:
                                        letter = seed_sequence[position] 
                                        mutations_seq += letter
                                        
                                pop[new_lin_name] = mutations_seq
                            else:
                                for position in range(0,len(seed_sequence)):
                                    if position in positions:
                                        letter = mutate(seed_sequence, position,qmatrix)
                                        mutations_seq += letter
                                    elif position in homo_mut:
                                        letter = seq_homo[position] 
                                        mutations_seq += letter
                                    else:
                                        letter = seed_sequence[position] 
                                        mutations_seq += letter
                        elif num_mut == 0:
                            for position in all_pos:
                                if position in homo_mut:
                                    letter = seq_homo[position]
                                    mutations_seq += letter
                                else:
                                    letter = seed_sequence[position] 
                                    mutations_seq += letter
                        pop[new_lin_name] = seed_sequence
                    else:
                        continue

                ## Recombination
                numRec = round(sum(rng.poisson(rec_rate, len(pop)))) 
                if numRec == 0: # Enable at least one inter-lineage recombinant and one intra-lineage recombination in each generation
                    numRec = 1
                already_id = list(pop.keys())
                # Inter-lineage recombination 
                for i in range(numRec):
                    flag_x1 = flag_x2 = 0
                    lin_name_x1 = lin_name_x2 = ""
                    # Select the sequences sourcing from different seeds(lineages), should exclude recombinant.
                    while flag_x1 == 0 or  flag_x2 == 0 or lin_name_x1 == lin_name_x2:
                        import random
                        x1_lin, x2_lin = random.sample(already_id, 2)
                        lin_name_x1,flag_x1 = get_lin_name(x1_lin)
                        lin_name_x2,flag_x2 = get_lin_name(x2_lin)

                    import random
                    pop = generate_recom_bk(seq_length, generation, x1_lin, x2_lin, pop)

                # Intra-lineage recombination 
                for i in range(numRec):
                    flag_xAA = 0
                    while flag_xAA == 0:
                        xAA1_lin = random.sample(already_id, 1)[0]
                        lin_name_xAA,flag_xAA = get_lin_name(xAA1_lin)
                    # find a sequence souced from the same seed(lineage)
                    random.shuffle(already_id)
                    for aa in already_id:
                        if ("bk" not in aa) and (lin_name_xAA in aa) and (aa != xAA1_lin):
                            xAA2_lin = aa
                            break
                    pop = generate_recom_bk(seq_length,generation,xAA1_lin,xAA2_lin,pop)
            ## The end of the simulation process and the simulted dataset in this independent trial has been generated.
            sampler_num = int(len(pop)/sample_rate)
            seq_samples_id = random.sample(pop.keys(), sampler_num)
            
            seq_samples = {}
            for key in pop:
                if key in seq_samples_id:
                    seq_samples[key] = pop[key]
            
            outpath = code_path+"gener"+str(genera)+"/turns"+str(turns)+"/"
            creat_dir(outpath)
            with open(outpath+"gener"+str(genera)+"_turns"+str(turns)+".fasta","w+") as h:
                for key in seq_samples:
                    h.write(">"+key+"\n")
                    h.write(seq_samples[key]+"\n")


# based on empirical n=1180 mutations counts from 7 countries
qmatrix = {
    'A': { "C": 0.1083, "G": 0.7000, "T": 0.1917 },
    'C': { "A": 0.0475, "G": 0.0033, "T": 0.9492 },
    'G': { "A": 0.2102, "C": 0.0931, "T": 0.6967 },
    'T': { "A": 0.1025, "C": 0.795, "G": 0.1025 }
    }

temp_file = code_path ="/home/soniali/Desktop/02_recom_230203/data/simulation_test/"
creat_dir(code_path)
seed_gen = 5
seed_gene_homo = 1 
num_seeds = 30 

turns_list = 5
sample_rate = 10 
parallel_prop = 0.1 
len_UAB = 4
seq_length = 29903
mut_rate = float(8E-4/6) # two months

seed_mut_perlin = 8 
generations_list = [6,7,8,9,10,11]  
rec_rate = 0.1

import random
simulator_withhomo(turns_list,seed_mut_perlin,generations_list, rec_rate)

