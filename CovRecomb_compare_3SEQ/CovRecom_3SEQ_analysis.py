import copy
import os
from numpy import *
import pandas as pd
import numpy as np
from Bio import SeqIO
import datetime 
import re
from scipy.stats import hypergeom
import copy
from statsmodels.sandbox.stats.multicomp import multipletests
import subprocess
os.chdir(os.getcwd()+"/")

def extract_ref_mut(ref_seq, seq_samples, lin_all_mutate, s):
    for g in range(len(seq_samples[s])):
        if seq_samples[s][g] != ref_seq[g]:
            mut = str(g+1)+"_"+str(seq_samples[s][g])
            lin_all_mutate.append(mut)
    return lin_all_mutate


def sample_lin_time(lin_A,lin_time):
    if (lin_A in lin_time):
        timeA = lin_time[lin_A]
    else:
        timeA = ''
    return timeA


def min_pairs(dic):
    if len(dic) == 0:
        return []
    min_val = min(map(lambda v: v[1], dic.items()))
    return [item for item in dic.items() if item[1] == min_val]


def sort_humanly(v_list): 

    def tryint(s):                       
        try:
            return int(s)
        except ValueError:
            return s

    def str2int(v_str):
        return [tryint(sub_str) for sub_str in re.split('([0-9]+)', v_str)]

    return sorted(v_list, key=str2int,reverse=False)


def calcul_bk_test(lin_A_draw,lin_B_draw,Lineage_v,epiV):
    feature_SNPA = Lineage_v[lin_A_draw]
    feature_SNPB = Lineage_v[lin_B_draw]

    A_B_shared = set(feature_SNPA) & set(feature_SNPB)
    UA_mutate = (set(feature_SNPA) & set(epiV)) - set(A_B_shared)
    UB_mutate = (set(feature_SNPB) & set(epiV)) - set(A_B_shared)

    UA_mutate_unique = []
    UB_mutate_unique = []
    
    lin_record = ""
    for j in epiV:
        if j in UA_mutate:
            UA_mutate_unique.append(j)
            lin_record = lin_record+"A"
        elif j in UB_mutate:
            UB_mutate_unique.append(j)
            lin_record = lin_record+"B"

    return lin_record,UA_mutate_unique,UB_mutate_unique


def unique_lin(aftertime_lin,Lineage_v,U_mutate_unique):
    num = 0
    same_ancient = []
    for can_A in aftertime_lin:
        if set(Lineage_v[can_A]) >= set(U_mutate_unique):
            same_ancient.append(can_A)
            num += 1

    if num ==1:
        return num
    elif num >1:
        same_A_num = mutation_lin_unique(same_ancient)
        if same_A_num == len(same_ancient):
            num =1
        else:
            num =2
    return num


def mutation_lin_unique(same_ancient):
    clean_same_ancient = []
    for mb in same_ancient:
        if "cluster" in mb:
            clean_same_ancient.append(mb.split("_")[1])
        else:
            clean_same_ancient.append(mb)

    same_num = 0
    for saan in clean_same_ancient:
        if len(saan.split("."))<3: 
            break
        else:
            same_01 = clean_same_ancient[0].split(".")[0:3]
            if saan.split(".")[0:3] != same_01:
                break
            elif saan.split(".")[0:3] == same_01:
                same_num +=1
    return same_num


def get_ch1_ch2_name(CH):
    CH1 = CH.split("bk")[0].split("_")[-1]
    if "_" in CH.split("bk")[1]:
        CH2 = CH.split("bk")[1].split("/")[0].split("_")[-1]
    else:
        CH2 = "lin"+CH.split("bk")[1].split("/")[0].split("lin")[-1]
    return CH1,CH2


def get_lin_name(epi):
    if "gen" in epi:
        CH1 = CH2 = "lin" + epi.split("_lin")[1]
    else:
        CH1 = CH2 = epi
    return CH1,CH2


def obtain_pattern(AB_epi, unique_A, unique_B):
    recom_pattern = ""
    for v in AB_epi:
        if v in unique_A:
            recom_pattern = recom_pattern + "X"
        elif v in unique_B:
            recom_pattern = recom_pattern + "Y"

    return recom_pattern

def bk_count(recom_pattern):
    start = recom_pattern[0]
    change = 0
    for R in recom_pattern:
        if R != start:
            change += 1
            start = R

    return change


def calcul_bk(lin_A_draw, lin_B_draw, Lineage_v, epiV):
    feature_SNPA = Lineage_v[lin_A_draw]
    feature_SNPB = Lineage_v[lin_B_draw]
    A_B_shared = set(feature_SNPA) & set(feature_SNPB)
    UA_mutate = (set(feature_SNPA) & set(epiV)) - set(A_B_shared)
    UB_mutate = (set(feature_SNPB) & set(epiV)) - set(A_B_shared)
    sample_special = set(epiV) - (set(feature_SNPA) | set(feature_SNPB))

    UA_mutate_unique = []
    UB_mutate_unique = []
    shared_mut = []
    denovo_mut = []

    lin_record = ""
    for j in epiV:
        if j in A_B_shared:
            shared_mut.append(j)
        elif j in UA_mutate:
            UA_mutate_unique.append(j)
            lin_record = lin_record + "X"
        elif j in UB_mutate:
            UB_mutate_unique.append(j)
            lin_record = lin_record + "Y"
        elif j in sample_special:
            denovo_mut.append(j)

    return lin_record, UA_mutate_unique, UB_mutate_unique, shared_mut, denovo_mut


def TimeConverter(ms):
    # 保留两位小数，但若ms太小，h就会显示为0。
    s = round(ms / 1000, 2)
    m = round(s / 60, 2)
    h = round(m / 60, 2)
    return s, m, h


def judge_yin(CH1,CH2,zhenyin,jiayin):
    if CH1 == CH2:
        zhenyin += 1
    else:
        jiayin += 1
    return  zhenyin,jiayin
            
def creat_dir(turns_file):
    import os
    if not os.path.exists(turns_file):
        os.makedirs(turns_file)
        

def recombination_detection_test(Strain_list_snp,len_UAB,linA_list,variants_all,feature_mutations,Lineage_v,mutaions_num, recom_total_num):
    linA_list_deep = copy.deepcopy(linA_list)
    zhenyin = 0
    jiayang = 0
    jiayin = 0
    zhenyang = 0
    correct_AB_num = 0
    false_AB_num = 0
    max_bk_num = 2
    
    for epi in Strain_list_snp:
        if "bk" in epi:
            CH1,CH2 = get_ch1_ch2_name(epi)
        else:
            CH1,CH2 = get_lin_name(epi)

        epi_record = {}
        epiV = variants_all[epi]
        epi_feat = len(set(epiV) & set(feature_mutations))

        ### P-value for Non-recombination
        for lin_A in linA_list_deep:
            all_AA = len(Lineage_v[lin_A])
            all_AA_epi = len(set(Lineage_v[lin_A]) & set(epiV))
            pVal = hypergeom.sf(all_AA_epi-1,mutaions_num,all_AA,epi_feat)
            epi_record[str(lin_A)+"_"+str(lin_A)] = pVal
        # the least p-value for the Non-recombinant
        min_AA = min(epi_record, key = epi_record.get)

        ### P-value for Recombinant (A+B/A+B+A)
        A_already = []
        for A in linA_list_deep:
            A_already.append(A)
            A_epi = set(Lineage_v[A]) & set(epiV)
            if len(A_epi) < len_UAB:
                continue
            else:
                afterA_linB = set(linA_list_deep) - set(A_already)
                for B in afterA_linB:
                    B_epi = set(Lineage_v[B]) & set(epiV)
                    
                    if len(B_epi) < len_UAB:
                        continue
                    else:
                        unique_A = A_epi - B_epi
                        unique_B = B_epi - A_epi
                        union_AB_set = set(Lineage_v[A]) ^ set(Lineage_v[B])
                        AB_epi = sort_humanly(list(union_AB_set & set(epiV)))
                        recom_pattern = obtain_pattern(AB_epi, unique_A, unique_B)
                        if len(recom_pattern) <= 1:
                            continue
                        else:
                            change = bk_count(recom_pattern)
                            if change > max_bk_num:
                                continue
                            else:
                                all_AB = len(set(Lineage_v[A]) | set(Lineage_v[B]))
                                all_AB_epi = len(set(set(Lineage_v[A]) | set(Lineage_v[B])) & set(epiV))
                                pVal = hypergeom.sf(all_AB_epi - 1, mutaions_num, all_AB, epi_feat)
                                epi_record[str(A) + "_" + str(B)] =  pVal

        raw_pvals = list(epi_record.values())
        rejected, p_adjusted, _, alpha_corrected = multipletests(raw_pvals, alpha=0.05, method='bonferroni', is_sorted=False, returnsorted=False)
    
        lin_adjP = {}
        for p in range(len(p_adjusted)):
            lin_pair = list(epi_record.keys())[p]
            lin_adjP[lin_pair] = p_adjusted[p]
            
        min_adjp_pair = min(lin_adjP, key = lin_adjP.get)
        if min_adjp_pair == min_AA or lin_adjP[min_adjp_pair] >= 0.05:
            zhenyin, jiayin = judge_yin(CH1, CH2, zhenyin, jiayin)
            # continue
        else:
            epiV = sort_humanly(epiV)
            lin_A_draw, lin_B_draw = min_adjp_pair.split("_")[0],min_adjp_pair.split("_")[1]
            lin_record, UA_mutate_unique, UB_mutate_unique, shared_mut, denovo_mut = calcul_bk(lin_A_draw, lin_B_draw, Lineage_v, epiV)
            if set(UA_mutate_unique + UB_mutate_unique) - set(Lineage_v[min_AA.split("_")[0]]) == set():
                zhenyin, jiayin = judge_yin(CH1, CH2, zhenyin, jiayin)
            else:
                if CH1 == CH2: # FP
                    jiayang += 1
                elif CH1 != CH2: # TP
                    zhenyang += 1
                    result_recom = {str(min_adjp_pair.split("_")[0]), str(min_adjp_pair.split("_")[1])}
                    if result_recom - {CH1, CH2} == set():
                        correct_AB_num += 1
                    else:
                        false_AB_num += 1

    if zhenyang+jiayin == recom_total_num:
        detection_rate = TPR = round(zhenyang / (zhenyang + jiayin),4)
        FPR = round(jiayang / (jiayang + zhenyin),4)
        FDR = round(jiayang/(jiayang+zhenyang),4)
        accurate_rate = round(correct_AB_num/(correct_AB_num+false_AB_num),4)
        return detection_rate, TPR, FPR, FDR,accurate_rate
    else:
        print("CovRecomb ERROR!")


turns_list = 5
generations_list = [6,7,8,9,10,11]
num_seeds = 30
len_UAB = 4
compare3SEQ = True
code_path ="/home/soniali/Desktop/02_recom_230203/data/simulation_test/"
os.chdir(code_path)

####### CovRecomb
with open(code_path+"EPI_ISL_402125.fasta", "r") as f:
    for record in SeqIO.parse(f, "fasta"):
        ref_seq = str(record.seq)


for genera in generations_list:
    recrd_results = {}
    for turns in range(1,turns_list+1):
        trial_name = "gener"+str(genera)+"_turns"+str(turns)
        outpath = code_path+"gener"+str(genera)+"/turns"+str(turns)+"/"
        fasta_file = outpath+trial_name+".fasta"
        recom_total_num = 0
        for seq_record in SeqIO.parse(fasta_file, "fasta"):
            if "bk" in seq_record.id:
                CH1,CH2 = get_ch1_ch2_name(seq_record.id)
                if CH1 != CH2:
                    recom_total_num += 1
        
        starttime = datetime.datetime.now()
        # read sequences
        seq_samples = {}
        for seq_record in SeqIO.parse(fasta_file, "fasta"):
            seq_samples[seq_record.id] = str(seq_record.seq)
        sample_size = len(seq_samples)
        # mutation calling
        Strain_list_snp = list(seq_samples.keys())
        variants_all = {}
        for s in seq_samples.keys():
            temp = []
            for g in range(len(seq_samples[s])):
                if seq_samples[s][g] != ref_seq[g]:
                    mut = str(g+1)+"_"+str(seq_samples[s][g])
                    temp.append(mut)
            variants_all[s] = temp
        
        # Extract feature mutations
        Lineage_v = {}
        for num in range(1,num_seeds+1):
            lin_all_mutate = [] # Record all the mutations within each lineage
            count_target_lin = 0 # Record the number of samples within each lineage, except recombinants
            for s in Strain_list_snp:
                if "bk" not in s and int(s.split("lin")[1]) == num:
                    count_target_lin += 1                  
                    lin_all_mutate =  extract_ref_mut(ref_seq, seq_samples, lin_all_mutate, s)
                elif "bk" in s and int(s.split("lin")[-1].split("/")[0]) == num:
                    count_target_lin += 1
                    lin_all_mutate =  extract_ref_mut(ref_seq, seq_samples, lin_all_mutate, s)
                
            fv75 = []
            for m in set(lin_all_mutate):
                if lin_all_mutate.count(m) >= (0.75*count_target_lin):
                    fv75.append(m)
            Lineage_v["lin"+str(num)] = fv75 

        # CovRecomb detection
        linA_list = list(Lineage_v.keys())
        feature_mutation = []
        for l in Lineage_v:
            for v in Lineage_v[l]:
                feature_mutation.append(v)
        
        feature_mutations = list(set(feature_mutation))
        mutaions_num = len(feature_mutations)
        detection_rate, TPR, FPR, FDR, accurate_rate = recombination_detection_test(Strain_list_snp,len_UAB,linA_list, variants_all,feature_mutations,Lineage_v,mutaions_num, recom_total_num)

        endtime = datetime.datetime.now()
        time_microsecond = (endtime - starttime).seconds * 1000 + (endtime - starttime).microseconds / 1000
        time_second, time_minute, time_hour = TimeConverter(time_microsecond)
        recrd_results[trial_name] = [detection_rate, TPR, FPR, FDR, accurate_rate, time_microsecond, time_second, time_minute, time_hour,sample_size]

    my_df = pd.DataFrame.from_dict(recrd_results, orient='index')
    my_df.columns = ["detection_rate","TPR","FPR","FDR","accurate_rate", "microseconds","seconds","minutes","hours","sample_size"]
    my_df.loc[ "gener"+str(genera)+"_mean"] = my_df.apply(lambda x: x.mean())
    my_df.to_csv(code_path+"gener"+str(genera)+"/"+"gener"+str(genera)+"_CovRecomb_bk.csv",index=True)


### 3seq
"./3seq -f mysequencefile.phy -d"
def get_parent_name(row_info,loc):
    loc_info = row_info[loc]
    if "bk" not in loc_info and "gen" not in loc_info:
        PA = loc_info
    elif "bk" not in loc_info and "gen" in loc_info:
        PA = loc_info.split("_")[-1]
    elif "bk" in loc_info:
        PA = "FALSE"
    return PA

for genera in generations_list:
    SEQ3_results = {}
    for turns in range(1,turns_list+1):
        trial_name = "gener"+str(genera)+"_turns"+str(turns)
        outpath = code_path+"gener"+str(genera)+"/turns"+str(turns)+"/"
        fasta_file = outpath+trial_name+".fasta"
        reference_file = code_path+"EPI_ISL_402125.fasta"
        
        recom_total_num = 0
        non_recom_total_num = 0
        sample_size = 0
        for seq_record in SeqIO.parse(fasta_file, "fasta"):
            sample_size += 1
            if "bk" in seq_record.id:
                CH1,CH2 = get_ch1_ch2_name(seq_record.id)
                if CH1 != CH2:
                    recom_total_num += 1
                else:
                    non_recom_total_num += 1
            else:
                non_recom_total_num += 1
        
        SEQ3_outpath = outpath + trial_name+"_3SEQ"
        creat_dir(SEQ3_outpath)
        starttime = datetime.datetime.now()
        import time

        subprocess.call (["cd /home/soniali/Downloads/3seq_build/ && \
                        echo y | ./3seq -f %s -id %s -t1 -d" % (fasta_file, SEQ3_outpath+"/"+trial_name)],shell=True)
    
        endtime = datetime.datetime.now()
        time_microsecond = (endtime -starttime).seconds * 1000 + (endtime -starttime).microseconds / 1000
        time_second, time_minute, time_hour = TimeConverter(time_microsecond)
        seq_samples = {}
        for seq_record in SeqIO.parse(fasta_file, "fasta"):
            seq_samples[seq_record.id] = str(seq_record.seq)
      
        for file in os.listdir(SEQ3_outpath+"/"):
            if file[-7:] == ".3s.rec":
                file_path = SEQ3_outpath+"/"+file
                with open(file_path,"r") as hseq:
                    rows = hseq.readlines()

        seq3_false_pos = 0
        seq3_true_pos = 0
        accurate_num = 0
        if len(rows) == 1: # 重组体为0
            overall_correct_rate = 0
        else:
            seq3_true_pos = 0
            for row in range(1,len(rows)):
                row_info  = rows[row].split("\t")
                CH = rows[row].split("\t")[2]
                pvalue = float(rows[row].split("\t")[9])
                if "bk" in CH and pvalue < 0.05:
                    CH1_r,CH2_r = get_ch1_ch2_name(CH)
                    if CH1_r != CH2_r:
                        seq3_true_pos += 1
                        PA = get_parent_name(row_info,0)
                        PB = get_parent_name(row_info,1)
                        if {CH1_r,CH2_r} - {PA,PB} == set():
                            accurate_num += 1
                        else:
                            continue
                    else:
                        continue
                elif "bk" not in CH and pvalue < 0.05:
                    seq3_false_pos += 1
                    
        detection_rate = TPR = round(seq3_true_pos/recom_total_num,4)
        FPR = round(seq3_false_pos / non_recom_total_num,4)
        FDR = round( seq3_false_pos/(seq3_false_pos+seq3_true_pos),4)
        accurate_rate = round(accurate_num/seq3_true_pos,4)
        SEQ3_results[trial_name] = [detection_rate, TPR, FPR, FDR, accurate_rate, time_microsecond, time_second, time_minute, time_hour,sample_size]

    my_3seq_df = pd.DataFrame.from_dict(SEQ3_results, orient='index')
    my_3seq_df.columns = ["detection_rate","TPR","FPR","FDR","accurate_rate","microseconds","seconds","minutes","hours","sample_size"]
    my_3seq_df.loc[ "gener"+str(genera)+"_mean"] = my_3seq_df.apply(lambda x: x.mean())
    my_3seq_df.to_csv(code_path+"gener"+str(genera)+"/"+"gener"+str(genera)+"_3SEQ.csv",index=True)


