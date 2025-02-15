#!/bin/python3.9.4
#Ian Outhwaite & Sukrit Singh, 2022
#032722
#run_all.py
#runs tests for MMR single-target protocol

import JS_singletarget_032722_exclude_mutants as JS
import pandas as pd
import numpy as np
import sys
import openpyxl

if __name__ == "__main__":
    
    ######################## USER DEFINED DATASET, should note correct number of rows to use (cols too if desired)

    # Note, dataset must have the same formatting as the excel file for the PKIS2 dataset.
    # User should add null values for some of the initial daata descriptions (ex: chemotype family) if not present in the dataset to ensure overal formatting is consistent.
    
    dataset_D=pd.read_excel("./2011_Davis_dataset_formatted_human_kinases_only.xlsx", engine='openpyxl', nrows=74)
                                
    print(dataset_D)

    ######################## Define names of targets NOT to include in the off-target profiles (ex: disease relevant mutants)

    Davis_names = ['ABL1(E255K)-phosphorylated',	'ABL1(F317I)-nonphosphorylated',	'ABL1(F317I)-phosphorylated',	'ABL1(F317L)-nonphosphorylated',	'ABL1(F317L)-phosphorylated',	'ABL1(H396P)-nonphosphorylated',	'ABL1(H396P)-phosphorylated',	'ABL1(M351T)-phosphorylated'	,'ABL1(Q252H)-nonphosphorylated',	'ABL1(Q252H)-phosphorylated',	'ABL1(T315I)-nonphosphorylated',	'ABL1(T315I)-phosphorylated',	'ABL1(Y253F)-phosphorylated',	'ABL1-nonphosphorylated',	'ABL1-phosphorylated' ,'BRAF(V600E)', 'EGFR(E746-A750del)'	,'EGFR(G719C)',	'EGFR(G719S)'	,'EGFR(L747-E749del, A750P)'	,'EGFR(L747-S752del, P753S)',	'EGFR(L747-T751del,Sins)',	'EGFR(L858R)'	,'EGFR(L858R,T790M)'	,'EGFR(L861Q)',	'EGFR(S752-I759del)',	'EGFR(T790M)','FGFR3(G697C)', 'FLT3(D835H)',	'FLT3(D835Y)'	,'FLT3(ITD)'	,'FLT3(K663Q)',	'FLT3(N841I)',	'FLT3(R834Q)', 'KIT(A829P)',	'KIT(D816H)',	'KIT(D816V)',	'KIT(L576P)',	'KIT(V559D)',	'KIT(V559D,T670I)',	'KIT(V559D,V654A)', 'LRRK2(G2019S)',
'MET(M1250T)',	'MET(Y1235D)', 'PIK3CA(C420R)',	'PIK3CA(E542K)',	'PIK3CA(E545A)'	,'PIK3CA(E545K)'	,'PIK3CA(H1047L)',	'PIK3CA(H1047Y)',	'PIK3CA(I800L)'	,'PIK3CA(M1043I)'	,'PIK3CA(Q546K)' ,'RET(M918T)',	'RET(V804L)',	'RET(V804M)']

    ######################## USER DEFINED SETTINGS
    
    prior_type = 2  #1=use a beta distribution, 2=use a poisson distribution
    bdist_alpha = 3 #The alpha value for the beta distribution
    bdist_beta = 1 #The beta value  for the beta distribution
    mu = 700 #The mu value for the poisson distribution. We find 200 and 700 to be good values to use; use both in different tests to compare selectivity profiles
    size = 100000 #The number of points to sample to generate the selectivity prior. We find 100,000 to be a decent number, it facilitates reproucibility without being too large.
    on_target_inhibition_threshold = 90 #the minimum on-target activity required for inhibitors. Try 90, 80, 70 as desired
    number_top_UW_scores_to_maximize = 5 #the top number_top_UW_scores_to_maximize unweighed combinations to select for the subsequent inhibitor concentration optimization step. We suggest using 3 or 5 - only in very rare instances would a combination that doesn't perform well in unweighed scoring perform exceptionally better in weighed scoring, and reducing the number of combinations to optimize speeds up processing time at this step linearly.
    noise_variance = 2.5 #the variance of the gaussian noise added to measurments
    max_combination_iter = 5 #the maximum number of inhibitor combinations to try (ex: 3 means using 3 inhibitors together in one mix)
    influence = 0.1 #the weight that is used to add an additional penalty to the highest off-target effects (in the [95, 100] bin). We advise using a value of 0.1.
    R1_step = 5
    R2_step = 1.1
    output_name = '' #Name of the output file
    num_p = 1 #number of proccesses to call during multiprocessing

    ####################### Default Prior Settings

    default_prior_settings = (prior_type, (bdist_alpha, bdist_beta, mu, size)) #default mu=700 Poisson

    ####################### RUN

    if len(sys.argv) > 1:
        print('Number of processes to use: ' + str(sys.argv[1]))
        num_p = int(sys.argv[1])

    ##################### FIRST DATASET
        
    dataset = dataset_D

    ########## Variable Prior Tests: Repeat each 5x 

    for x in range(0,5):
    
        excel_prefix = '032822_Davis_90thresh_x'+str(x)+'_'

        for new_mu in [200,700]:
            output_name = excel_prefix + 'poisson_' + str(new_mu)
            new_prior_type = 2
            new_prior_settings = (new_prior_type, (bdist_alpha, bdist_beta, new_mu, size))
            JS.main_function(dataset, new_prior_settings, noise_variance, max_combination_iter, influence, output_name, num_p, on_target_inhibition_threshold, number_top_UW_scores_to_maximize, Davis_names, R1_step, R2_step)


    
