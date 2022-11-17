
### Genetic variants associated with long-terminal repeats can diagnostically classify Cannabis varieties (Oultram et al., 2022).             ###
### This script is used to model the fragmentation of genomic DNA in the form of a fasta file based on specified restriction enzyme sequences ###
 

### Input fasta file location and output file destination ###
fastafile_input = ''
fragmentfile_output = ''


### Determine upper and lower thresholds for fragments to print to file ###
fragment_print_max_limit = 600 # maximum size of fragments to print out
fragment_print_min_cutoff = 50  # min size of fragments to print out


### Determine upper and lower thresholds for fragments to graph ###
graph_upper_limit = 600 # maximum size of fragments to print out
graph_lower_limit = 50  # min size of fragments to print out


### Create a dictionary of restriction enzyme sequences to digest your fasta file ###

re_dictionary = {#"HhaI":"GCGC", 
                 #"BstNI":"CC(A|T)GG",  
                 #"EcoRI":"GAATTC",
                 #"PstI":"CTGCAG",
                 #"AluI":"AGCT"
                 #"HindIII":"AAGTCC",
                 #"SacI":"GAGCTC"
                   }


output_file_object = open(fragmentfile_output, 'w')


### Convert RE dictionary to list of RE sites ###
re_list = []
for eachentry in re_dictionary:
    re_list.append(re_dictionary[eachentry])
                

with open(fastafile_input , 'r') as fastafile:

    ### Dictionary to hold RE fragment sizes and frequency ###
    datadictionary = {}  
    

  
    for line1 in fastafile:   
                             
        ### Check for fasta name in first line ###
        if line1.startswith(">"):
            fastaheader = line1[1:].rstrip() #### Drop '>' and remove newline ###
            olddna = ''
            dna_position_start = 0  ### Variable to hold the base position in the fasta entry we've processed upto, for printing fragments ###
            print(fastaheader)
            continue
            
           
        concatDNA = olddna + line1.rstrip().upper()  ### Pull out next line of DNA, append to previous DNA for which there was no RE site match ###                                    


        re_results = []
        for each_sequence in re_list:
            if concatDNA.find(each_sequence) != -1:  ### If no RE match, then skip to next sequence ###
                re_results.append(concatDNA.find(each_sequence))  ### If RE site Found, record occurance in list ###
            
        ### If no RE site found, copy current DNA to 'olddna' variable for next loop ###
        if not re_results:  ### Check if list is empty ###
            olddna = concatDNA  

        ### If RE site is found, print out fragment and update data dictionary ###
        else:
            ### Add 2 to size and slice, or else will pick up same RE site in next loop ###
            fragmentsize = min(re_results) + 2
            
            ### Info for printing to fasta file ###
            fragment_start = dna_position_start
            fragment_end = dna_position_start + fragmentsize
            
            re_fragment = concatDNA[:fragmentsize] ### DNA sequence of RE fragment ###
            
            dna_position_start += fragmentsize  ### Update variable with new position start coordinate ###
            olddna = concatDNA[fragmentsize:]   ### Read DNA that wasn't cut into old DNA variable for next loop ###

            ###############################################
            ###############################################
            ## Update data dictionary with fragment info ###
            if graph_upper_limit >= fragmentsize >= graph_lower_limit :
        
                    if fragmentsize in datadictionary:
                        datadictionary[fragmentsize] += 1
                    else:
                        datadictionary[fragmentsize] = 1

            ###############################################
            ###############################################    
            ### If fragment is less than or equal to limit size, print fragment to file for downstream mapping ###
            if fragment_print_max_limit >= fragmentsize >= fragment_print_min_cutoff:
                
     
        
                ### Print fragments to file ###
                name1 = " ".join( [">", str(fragment_start), str(fragment_end), str(fragmentsize), fastaheader]   )
                output_file_object.write(name1 + "\n")
                output_file_object.write(re_fragment  + "\n")
            



output_file_object.close()



#%%  line plot of data dictionary

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy.signal import savgol_filter

### Convert datadictionary to pandas df ###
df1 = pd.DataFrame(datadictionary, index=[0])  
a2 = df1.iloc[0,]

### Extract out entries into lists ###
sizex = []
massy = []

for eachentry in datadictionary:
    sizex.append(eachentry)
    massy.append(datadictionary[eachentry]*eachentry)


### Savitzky-Golay filter to smooth graph lines to replicate LabChip output ###
#y_fil = savgol_filter(massy, 71, 5)
# massdf = pd.DataFrame()
# massdf['Size'] = sizex
# massdf['Mass'] = massy
# massdf.to_csv('', index=False, header=False)

### Can change the y axis to y=massgraph using below to change the y axis from base pair count per fragment size to count of fragment size ###

# massgraph = [m/s for m, s in zip(map(int, massy), map(int, sizex))]
# print(massgraph)

#g = sns.lineplot(x=sizex, y=y_fil, size=, color='', legend=False)

### Graph features that can be altered as needed ###
#plt.xlabel('')
#plt.ylabel('')
#plt.tight_layout()
#plt.xlim(, )
#plt.ylim(, )



#plt.savefig('...png', dpi=)
#plt.show()

### Save the dataframe as a csv then used $ sed 1d *.csv > merged.csv to merge them all ###
#a2.to_csv('/Users/...csv', index=True)

#%%  line plot of df
import pandas as pd
import seaborn as sns
import matplotlib as plt

### Read a csv generated from the merged output from above. The csv is handled in excel prior to this to make it readable given the column headers etc ###
csv = pd.read_csv(r'')


### Alter graph features as needed ###
#sns.set_palette("")
#sns.set_style("", {"axes.facecolor": ".9"})
#sns.despine()
#sns.set_context('', font_scale=)
### x and y values relate to previous csv file amended in excel or similar ###
#res = sns.regplot(x='', y='',
                  #data=csv,
                  #lowess=True,
                  #scatter=False,
                  #linewidth=0.7,
                  #ci='sd',)
                  #kind='line',
                  #size=0.1,
                  #height=6,
                  #aspect=2)
#(res.set_axis_labels("Size [bp]")
#   .tight_layout(w_pad=0))
#res.set(xscale='') 
#res.fig.suptitle("", y=)



#res.set(ylim=(, ))
#res.set(xlim=(, ))
#res.savefig('', bbox_inches="tight", dpi=1500)
#plt.show()
