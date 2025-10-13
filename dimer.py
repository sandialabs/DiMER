import os
import sys
import shutil
import subprocess
import re
import time
import pickle
#also pandas but the import is down lower
from datetime import date
start = time.time()


inputs = sys.argv

#basic example:
#python main.py -in <input_path> -out <output path>

print("""
 _____      __  __  _____  _____ 
|  __ \\  _ |  \\/  ||  ___||  __ \\ 
| |  | |(_)| \\  / || |__  | |__) |
| |  | | _ | |\\/| ||  __| |  _  / 
| |__| || || |  | || |___ | | \\ \\ 
|_____/ |_||_|  |_||_____||_|  \\_|
""",flush=True)

print('\nVersion 1.0', flush=True)
print("Today's Date: " + str(date.today()), flush=True)


if "-h" in inputs:
    print("\n\n\n#############################################################################")
    print("DiMER")
    print("Created Dec 2024. Last Updated 06/18/2025.")
    print("#############################################################################\n")
    
    # Tool description
    print("This tool serves to improve the ease of querying and comparing a set of")
    print("predicted translated ORFs against several pertinent HMM and Blast databases.")
    print("You can add additional databases by modifying `diamond_db_paths.txt` or `hmmer3_db_paths.txt`")
    print("to include the path to new databases.")
    print("Please feel free to reach out to eltorra@sandia.gov for questions and assistance.\n")
    
    print("############################################################################################")
    print("############################################################################################")
    
    #####################
    # Required parameters
    print("Required:")
    print("############################################################################################\n")
    print("\t-in  <specify path to file of translated predicted coding sequences in fasta format>\n")
    print("\t-out <specify path to output directory>\n")
    print("############################################################################################\n")
    
    #####################
    # Optional parameters
    print("Optional Parameters:")
    print("############################################################################################\n")

    print("   -clean")
    print("        Deletes all created folders and files, except for the final output `.tsv`.\n")

    print("   -restart")
    print("        Restarts the run. The script will not regenerate HMM and BLAST outputs that already exist")
    print("        in the output directory unless this flag is used.\n")
    print("        **WARNING**: This flag will remove everything in the specified output directory.")
    print("        Please use a unique directory for this program.\n")

    print("   -unique_name <new output file name>")
    print("        Specify a custom file name instead of `combined_output.tsv`.\n")

    print("   -best_hits")
    print("        Outputs a file called `best_hits.txt` (or `<unique_name>.best_hits.txt` if a unique name is specified).")
    print("        This file contains annotations with the highest BLAST score, prioritizing high-scoring hits")
    print("        (above the e-value threshold) that are not labeled as placeholder annotations.\n")

    print("   -replace <specify 'force' or 'gentle'>")
    print("        Replace placeholder annotations in `.prot` files with new annotations.")
    print("        - **Gentle**: Only replaces annotations if the strings from placeholders.txt appears in the line and can be replaced by a non-placeholder option with a sufficiently low e-value")
    print("        - **Force**: Reannotates all lines with the best scoring non-placeholder hit across all models (unless a non-placeholder option is unavailable then it will default to the best placeholder option)\n")
    print("        Both options create a new file with your original file prefix followed by `.updated.faa`")
    print("        in the specified output directory.\n")

    print("   -eval <e-value threshold>")
    print("        Specify the maximum e-value threshold for inclusion in the list of best hits for final annotation")
    print("        of `.prot` files when using the `-replace` flag.")
    print("        - **Effect**: Using this flag without `-replace` has no effect.")
    print("        - **Default**: `1e-10`.\n")

    print("   -cpu <number of CPUs>")
    print("        Specify the number of additional CPUs to use for HMMER.")
    print("        - **Default**: 0 (no additional CPUs).")
    print("        - **Requirement**: Must be an integer.\n")

    print("   -blast_paths <path to file>")
    print("        Specify the path to a file containing paths to Diamond-formatted BLAST databases.")
    print("        - **Use Case**: Useful for specifying different sets of BLAST databases for different taxa in analysis pipelines.")
    print("        - **Default**: `diamond_db_paths.txt`.\n")

    print("   -hmm_paths <path to file>")
    print("        Specify the path to a file containing paths to HMMER3-formatted HMM databases.")
    print("        - **Use Case**: Useful for specifying different sets of HMM databases for different taxa in analysis pipelines.")
    print("        - **Default**: `hmmer_db_paths.txt`.\n")

    print("   -placeholder <path to file>")
    print("        Specify the path to a plain text file containing a list of placeholder annotations.")
    print("        - **Use Case**: Useful for specifying different placeholder annotations for different taxa.")
    print("        - **Default**: `placeholders.txt`.\n")

    ################
    # Additional instructions
    print("############################################################################################")
    print("If you want to add new HMM or Diamond-formatted BLAST datasets: run `python repickle.py`")
    print("after you've altered the `Databases_key.tsv` file to include new descriptions for HMM models you've")
    print("added. You will only need to do this once after altering this file. You will also need to")
    print("add the appropriate paths to your databases in the files `hmmer3_db_paths.txt` (HMM formatted")
    print("with HMMER3) and `diamond_db_paths.txt` (BLAST paths formatted with Diamond). Make sure the")
    print("descriptors for the paths in these files match with those in the first column of `Databases_key.tsv`.")
    print("############################################################################################\n")

    print("(I hope you have a nice day)\n\n")
    
    quit()

#check user has required software
#############################################
#check pandas is installed and initialize
#############################################
try:
    import pandas as pd
except ImportError:
    print("\nThe Pandas Python library is not installed or found. Please make sure it is installed.\n", flush=True)
    exit()
#############################################
#Check Diamond is Installed
#############################################
def check_diamond_installed():
    try:
        result = subprocess.run(['diamond', '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if result.returncode == 0:
            print("\n###########\n\n\nDiamond is installed.", flush=True)
            print("Version information:", flush=True)
            print(result.stdout, flush=True)
            print('\n###########\n', flush=True)
        else:
            print("Diamond is not installed or is not found. Please install Diamond to use this program. https://github.com/bbuchfink/diamond", flush=True)
            print(result.stderr, flush=True)
            exit()
    except FileNotFoundError:
        print("Diamond is not installed or is not found. Please install Diamond to use this program. https://github.com/bbuchfink/diamond", flush=True)
        exit()
check_diamond_installed()
#############################################
#Check HMMER3 is Installed
#############################################
def check_hmmer_installed():
    try:
        result = subprocess.run(['hmmsearch', '-h'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if result.returncode == 0:
            print("\n###########\n\n\nHMMER3 is installed.", flush=True)
            print("Version information:", flush=True)
            stuff = result.stdout.splitlines()
            print(stuff[1].strip('# ') + '\n', flush=True)
            print('\n###########\n', flush=True)
        else:
            print("HMMER3 is not installed or is not found. Please install HMMER3 to use this program. https://github.com/EddyRivasLab/hmmer", flush=True)
            print(result.stderr, flush=True)
            exit()
    except FileNotFoundError:
        print("HMMER3 is not installed or is not found. Please install HMMER3 to use this program. https://github.com/EddyRivasLab/hmmer", flush=True)
        exit()
check_hmmer_installed()
###############################################
#handles input
###############################################
cds_path = ''
if "-in" in inputs:
    i = inputs.index("-in")
    cds_path = inputs[i+1]
    #if cds_path[-1] == "/":
    #    pass
    #else:
    #    cds_path += "/"
    print("Input CDS (.faa) file: " + cds_path, flush=True)
elif not "-in" in inputs: #exit if no input given
    print("FATAL ERROR: User must specify path to predicted coding sequences in '.faa' format with the flag '-in <path-to-file>'", flush=True)
    quit()

###############################################
#Check if CDS file is translated
###############################################
def is_dna_sequence(sequence):
    dna_bases = set("ATCGN")
    return all(char in dna_bases for char in sequence)

def check_first_fasta_sequence(fasta_filename):
    with open(fasta_filename, 'r') as fasta_file:
        sequence = ""
        for line in fasta_file:
            line = line.strip()
            if line.startswith(">"):
                if sequence:
                    # We have already read the first sequence, so we stop here
                    break
            else:
                sequence += line
        # Check the first sequence
        if sequence and is_dna_sequence(sequence):
            return True
    return False

# Example usage
fasta_file = cds_path
if check_first_fasta_sequence(fasta_file):
    print("The first sequence in the FASTA file is a DNA sequence. All sequences must be translated. Exiting program.")
    quit()



###############################################
#handles output and restart
###############################################
out_dir_path = ''
restart = False
if "-out" in inputs:
    i = inputs.index("-out")
    out_dir_path = inputs[i+1]
    if out_dir_path[-1] == "/":
        pass
    else:
        out_dir_path += "/"
    if os.path.exists(out_dir_path):
        if "-restart" in inputs:
            try:
                shutil.rmtree(out_dir_path)
                os.mkdir(out_dir_path)
                restart=True
            except OSError as e:
                print("Error, unable to restart: %s - %s." % (e.filename, e.strerror), flush=True)
                print('Exiting. You may try restarting after manually deleting this folder or specifying a new output directory.', flush=True)
                quit() # may be issue with permissions of sym links
    else:
        try:
            os.mkdir(out_dir_path)
        except Exception as e:
            print("Exiting. Unable to create output directory due to the following system error: ", flush=True)
            print(str(e), flush=True)
            quit()
elif not "-out" in inputs:
    out_dir_path += './dimer_out/'
    try:
        os.mkdir(out_dir_path)
    except Exception as e:
        print('Exiting. Unable to write to current directory due to the following system error: ', flush=True)
        print(str(e), flush=True)
        quit()
print("Annotater output directory: " + out_dir_path)

###############################################
#Handles cpu
###############################################
cpu = 0
if '-cpu' in inputs:
    i = inputs.index("-cpu")
    #print(inputs[i+1])
    try:
        #print(inputs[i+1])
        cpu += int(inputs[i+1])
    except Exception as e:
        print('-cpu needs to be an integer, please.', flush=True)
        #print(e)
        quit()
print("Number of CPUs specified (note, too many can actually slow the program down!): " + str(cpu), flush=True)
###############################################
#handles clean
###############################################
clean = 'TRUE' if "-clean" in inputs else 'FALSE'

###############################################
#handles combined_output unique name
###############################################
unique_name = 'combined_output.tsv'
best_hits = 'best_hits.txt'
if "-unique_name" in inputs:
    i = inputs.index("-unique_name")
    unique_name = inputs[i+1]
    if not '.' in unique_name:
        unique_name = unique_name + '.tsv'
    if "-best_hits" in inputs:
        best_hits = unique_name + '_' + best_hits

output_file = os.path.join(out_dir_path, unique_name)
print('Final output will be written to: ' + str(output_file), flush=True)
################################################
# parse replace options
################################################
replace = 'FALSE'
if "-replace" in inputs:
    i = inputs.index("-replace")
    #replace = 'TRUE'
    replace = inputs[i+1]
    if replace == "force":
        print('-replace = force. All annotations will be replaced with new annotations.', flush=True)
        pass
    elif replace == "gentle":
        print('-replace = gentle. All known placeholder annotations will be replaced with new annotations.', flush=True)
        pass
    else:
        print('Sorry to interrupt, but you used the flag "-replace" without specifying "force" or "gentle". This is case sensitive!', flush=True)
        quit()

################################################
# handles eval threshhold
################################################
e_value_threshold = 1e-10
if "-eval" in inputs:
    i = inputs.index("-eval")
    e_value_threshold = inputs[i+1]
    try:
        float(e_value_threshold)
    except:
        print('Sorry to interrupt, but you used the flag "-eval" without giving a float!', flush=True)
        quit()

#######################################
#gets locaion of dimer.py
#######################################
script_path = os.path.abspath(__file__)
dimer_directory = os.path.dirname(script_path)

################################################
# handles alternate blast or hmmer text files
################################################
blast_paths = dimer_directory + '/diamond_db_paths.txt'
hmm_paths = dimer_directory + '/hmmer3_db_paths.txt'
if "-blast_paths" in inputs:
    i = inputs.index("-blast_paths")
    blast_paths = inputs[i+1]
if "-hmm_paths" in inputs:
    i = inputs.index("-blast_paths")
    hmm_paths = inputs[i+1]

################################################
# handles hypothetical synonyms file and alternate hypothetical synonyms text file
################################################
placeholder = dimer_directory + '/placeholders.txt'
if "-placeholder" in inputs:
    i = inputs.index("-placeholder")
    placeholder = inputs[i+1]
file = open(placeholder, 'r')
hyp_syn_list = []
for line in file:
    if not line.startswith('#'):
        hyp_syn_list.append(line.strip('\n'))
file.close()
hyp_syn_str= '|'.join(hyp_syn_list)



#################################
#Run HMMER3
#################################
hmmer_todo = {}
file = open(hmm_paths, 'r')
for line in file:
    if not line.startswith('#') and not line.startswith('\n'):
        try:
            name = line.strip(' ').strip('\n').split('=') [0]
            hmmer_todo[name] = ''
            path = line.strip(' ').strip('\n').split('=') [1]
            if path.startswith('.'):
                path = path.replace('.', dimer_directory, 1)
            hmmer_todo[name] += path
        except:
            pass

file.close()

print('\nRunning HMMER (This is not a drill)...', flush=True)
for name in hmmer_todo:
    hmm_out_path = out_dir_path + name + '/'
    database = str(hmmer_todo[name])
    print('Comparing Query to HMM database: ' + database, flush=True)
    result = out_dir_path + name + '/' + name + '.out'
    result2 = out_dir_path + name + '/' + name + '.tsv'
    if not os.path.exists(hmm_out_path):
        os.mkdir(hmm_out_path)
        print('Comparing Query to HMMER database: ' + database, flush=True)
        if not os.path.isfile(result):
            os.system('hmmsearch ' +'--cpu '+str(cpu)+' --tblout ' + result2 + ' ' + database + ' ' + cds_path +  ' > '+result)
    else:
        if not os.path.isfile(result):
            os.system('hmmsearch '+'--cpu '+str(cpu)+ '--tblout ' + result2  + ' ' + database + ' ' + cds_path + ' > '+result)
        elif os.stat(result).st_size==0:
            if restart is True:
                os.remove(result)
                os.system('hmmsearch ' +'--cpu '+str(cpu)+'--tblout ' + result2 + ' ' + database + ' ' + cds_path +  ' > '+result)
        else:
            pass

#################################
# Run Diamond on BLAST DATABASES
#################################
blast_todo = {}
file = open(blast_paths, 'r')
for line in file:
    if not line.startswith('#') and not line.startswith('\n'):
        try:
            name = line.strip(' ').strip('\n').split('=') [0]
            blast_todo[name] = ''
            path = line.strip(' ').strip('\n').split('=') [1]
            if path.startswith('.'):
                path = path.replace('.', dimer_directory, 1)
            blast_todo[name] += path
        except:
            pass
file.close()

print('\nSo anyway, I started BLASTing...', flush=True)
for name in blast_todo:
    blast_out_path = out_dir_path + name + '/'
    database = str(blast_todo[name])
    result = out_dir_path + name + '/' + name + '.tsv'
    print('Comparing Query to BLAST database: ' + database, flush=True)
    if not os.path.exists(blast_out_path):
        os.mkdir(blast_out_path)
        if not os.path.isfile(result):
            os.system('diamond blastp -d ' + database + ' -q ' + cds_path + ' -o ' + result)
    else:
        if not os.path.isfile(result):
            os.system('diamond blastp -d ' + database + ' -q ' + cds_path + ' -o ' + result)
        elif os.stat(result).st_size==0:
            if restart is True:
                os.remove(result)
                os.system('diamond blastp -d ' + database + ' -q ' + cds_path + ' -o ' + result)
        else:
            pass

#################################
#Parse Input File
#################################
input_dict = {}

with open(cds_path, 'r') as file:
    for line in file:
        if line.startswith('>'):
            try:
                target_name = ((line.split() [0]).replace('>', '')).strip('\n')
                annotation = (line.split(target_name) [1]).strip('\n')
                match = re.search(r'\[protein=([^\]]+)\]', annotation) #extract annotation from complex NCBI headers... doesn't retain the header though...should we?
                if match:
                    annotation = match.group(1)
                input_dict[target_name] = annotation
                
            except:
                target_name = ((line.split() [0]).replace('>', '')).strip('\n')
                annotation = 'NA'
                input_dict[target_name] = annotation
                


############################################
#PARSE HMM OUTPUT
############################################
def parse_hmm_data(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    target_data = {}

    for line in lines:
        if not line.startswith('#'):
            columns = re.split(r'\s+', line.strip()) # this doesn't always work because when names are long, sometimes re ###is only one whitespace
            target_name = columns[0]
            query_name = columns[2]

            e_value = float(columns[4])

            columns = re.split(r'\s{2,}', line.strip()) # grabs last column of split on 2+ whitespace
            description = columns[-1]

            if target_name not in target_data or e_value < target_data[target_name]['e_value']:
                target_data[target_name] = {
                    'e_value': e_value,
                    'description': description,
                    'query_name': query_name,
                    }

    return target_data
###############################################

hmmer_final={}
for name in hmmer_todo:
    print('Parsing HMMER output: ' +name, flush=True)
    hmmer_final[name] = parse_hmm_data(out_dir_path + name + '/' + name + '.tsv')
#print(hmmer_final)


############################################
#PARSE BLAST OUTPUT
############################################
def parse_blast_data(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    target_data = {}

    for line in lines:
        if not line.startswith('#'):
            columns = line.strip('\n').split('\t')
            target_name = columns[0]
            e_value = float(columns[10])
            description = columns[1]

            if target_name not in target_data or e_value < target_data[target_name]['e_value']:
                target_data[target_name] = {
                    'e_value': e_value,
                    'description': 'NA',
                    'query_name': description,
                    }

    return target_data
###############################################

blast_final = {}
for name in blast_todo:
    print('Parsing Diamond output: ' +name, flush=True)
    blast_final[name] = parse_blast_data(out_dir_path + name + '/' + name + '.tsv')
#print(blast_final)



#################################################
#FINAL OUTPUT
###########################################
# Collect all unique target names
# Collect all unique target names
all_target_names = set(input_dict.keys())
for d in hmmer_final.values():
    all_target_names.update(d.keys())
for d in blast_final.values():
    all_target_names.update(d.keys())

with open(dimer_directory + '/Databases_key.pickle', 'rb') as handle:
    b = pickle.load(handle)

# Write to TSV file
#making header
with open(output_file, 'w') as f:
    header = ['Target Name', 'Input Annotation'] 
    for db_name in hmmer_todo.keys():
        header.extend([f'HMMER {db_name} E-value', f'HMMER {db_name} Query_Name', f'HMMER {db_name} Main Description', f'HMMER {db_name} Generalized Description'])
    for db_name in blast_todo.keys():
        header.extend([f'BLAST {db_name} E-value', f'BLAST {db_name} Query_Name', f'BLAST {db_name} Main Description',f'BLAST {db_name} Generalized Description'])
    f.write('\t'.join(header) + '\n')

    # Write data
    for target_name in sorted(all_target_names):
        row = [target_name, input_dict.get(target_name, 'NA')]
        for db_name in hmmer_todo.keys():
            hmmer_data = hmmer_final.get(db_name, {}).get(target_name, {'e_value': 'NA', 'query_name': 'NA', 'descriptor': 'NA', 'gen_descriptor': 'NA'})
            #print(b[db_name])
            try:
                row.extend([hmmer_data['e_value'], hmmer_data['query_name'], b[db_name][(hmmer_data['query_name'])]])
            except:
                row.extend([hmmer_data['e_value'], hmmer_data['query_name'], 'NA', 'NA'])
        for db_name in blast_todo.keys():
            blast_data = blast_final.get(db_name, {}).get(target_name, {'e_value': 'NA', 'query_name': 'NA', 'description': 'NA', 'gen_description': 'NA'})

            try:
                row.extend([blast_data['e_value'], blast_data['query_name'], b[db_name][(blast_data['query_name'])]])
            except:
                #row.extend([hmmer_data['e_value'], hmmer_data['query_name'], 'NA', 'NA'])
                row.extend([blast_data['e_value'], blast_data['query_name'], 'NA', 'NA'])
        f.write('\t'.join(map(str, row)) + '\n')

print(f"Final output written to {output_file}", flush=True)
print("\n###################################\n", flush=True)
#################################################
#function for checking that no description contains values from "hypothetical synonyms list"
def is_valid_description(description):
    if pd.isna(description):
        return False
    
    description_lower = description.lower()
    
    for term in hyp_syn_list:
        term_lower = term.lower().strip()
        if '#' in term_lower: #replaces # with any string of integers
            pattern = term_lower.replace('#', r'\d+')
            if re.search(pattern, description_lower):
                return False

        if term_lower in description_lower:
            return False
            
    return True
#################################################
#Update prot file GENTLE (only replace if annotation is hypothetical/unannotated/unknown)
#################################################
if replace == 'gentle':
    df = pd.read_csv(output_file, sep='\t')

    # Filter rows where column 2 contains 'hypothetical protein' (case insensitive) etc

    df.iloc[:, 1] = df.iloc[:, 1].fillna('').astype(str)
    filtered_df = df[df.iloc[:, 1].str.contains(hyp_syn_str, case=False) | df.iloc[:, 1].isin(['NA', '', 'NaN', 'nan'])] 
    total_hyp_entry = (len(filtered_df))

    results = []
    # Iterate through the filtered rows
    for index, row in filtered_df.iterrows():
        min_e_value = float('inf')
        corresponding_description = None

        # Check for E-value columns
        for col in df.columns:
            if 'E-value' in col:
                e_value = row[col]
                try:
                    e_value = float(e_value)  
                except (ValueError, TypeError):
                    continue  

                # Determine the corresponding description
                if col.startswith("HMMER"):
                    description = row[df.columns[df.columns.get_loc(col) + 2]]
                    if description == 'NA': #if there is no additional info provided by the annotation tsv, just retrieve the HMM or BLAST ID
                        description = row[df.columns[df.columns.get_loc(col) + 1]]
                elif col.startswith("BLAST"):
                    description = row[df.columns[df.columns.get_loc(col) + 2]]
                    if description == 'NA': #if there is no additional info provided by the annotation tsv, just retrieve the HMM or BLAST ID
                        description = row[df.columns[df.columns.get_loc(col) + 1]]

                # Check if the description does not contain 'hypothetical protein', 'unknown function', or 'unannotated protein' etc
                if not is_valid_description(description) == False:
                        # Update minimum E-value and corresponding description if current E-value is lower
                    if e_value < min_e_value:
                        min_e_value = e_value
                        corresponding_description = description
        
        # Store the result for the current row if a valid E-value was found
        if min_e_value < e_value_threshold:
            results.append({
                'Row Index': index,
                'First Column Value': row[df.columns[0]], 
                'Minimum E-value': min_e_value,
                'Description': corresponding_description
            })
    #print(results)
    if results:
        print("Gentle mode replaced annotations for: \n")
        for result in results:
            print(f"Row Index: {result['Row Index']}, First Column Value: {result['First Column Value']}, Minimum E-value: {result  ['Minimum E-value']}, Description: {result['Description']}\n", flush=True)
    else:
        print('No matching rows found with E-value < ' + str(e_value_threshold) + ' thus, no replacements were made.', flush=True)
    #calc # replaced
    hyp_replaced = len(results)
    print('Gentle mode improvement stats: ', flush=True)
    if not total_hyp_entry == 0:
        improvement = round(((hyp_replaced/total_hyp_entry)*100), 1)
        print('Total ' + '"' + '", "'.join(term for term in hyp_syn_list) + '" :' + str(total_hyp_entry), flush=True)
        print('Total replaced: ' + str(hyp_replaced), flush=True)
        print('Percent Improvment (%): ' + str(improvement), flush=True)
        print("\n###################################\n", flush=True)

    # Path to the FASTA file
    fasta_file_path = cds_path 

    base_name = os.path.splitext(os.path.basename(fasta_file_path))[0] 
    new_fasta_file = os.path.join(out_dir_path, f"{base_name}.updated.faa") 

    df = pd.read_csv(output_file, sep='\t')
    
    with open(fasta_file_path, 'r') as fasta_file:
        fasta_lines = fasta_file.readlines()

    updated_fasta_lines = []
    for line in fasta_lines:
        if line.startswith('>'):  
            for result in results:
                first_column_value = result['First Column Value']
                if first_column_value in line:  
                    line = '>' + first_column_value + ' ' + result['Description'] + '\n'
                    break  
        updated_fasta_lines.append(line)

    with open(new_fasta_file, 'w') as new_file:
        new_file.writelines(updated_fasta_lines)

    print("\n###################################\n", flush=True)
    print(f"Updated FASTA file written to: {new_fasta_file}", flush=True)
    print("\n###################################\n", flush=True)


#################################################
#Update prot file FORCE (replaces all)
###########################################

if replace == 'force' or '-best_hits' in inputs:
    besthit_output = os.path.join(out_dir_path, best_hits)
    df = pd.read_csv(output_file, sep='\t')
    try:
        # Filter out rows based on the hyp_syn_str or where the second column is 'NA', '', or 'NaN'
        df.iloc[:, 1] = df.iloc[:, 1].fillna('').astype(str)
        filtered_df_OG = df[df.iloc[:, 1].str.contains(hyp_syn_str, case=False, na=False) | df.iloc[:, 1].isin(['NA', '', 'NaN', 'nan'])]
        total_hyp_entry = len(filtered_df_OG)
    except:
        total_hyp_entry = len(df)
    
    results = []
    count = 0
    all_gene=0

    for index, row in df.iterrows():
        min_e_value = float('inf')
        corresponding_description = None
        corresponding_colid = None
        
        # Variables to track the best fallback option if no valid descriptions are found
        fallback_min_e_value = float('inf')
        fallback_description = None
        fallback_colid = None

        for col in df.columns:
            if 'E-value' in col:
                e_value = row[col]
                try:
                    e_value = float(e_value)
                except (ValueError, TypeError):
                    continue  # Skip to the next column if the E-value can't be converted to float

                # Determine the corresponding description
                description_col_index = df.columns.get_loc(col) + 2
                description = row[df.columns[description_col_index]]
                if description == 'NA':  # Fallback to the ID if no description
                    description = row[df.columns[description_col_index - 1]]

                # Check if the description is valid and the E-value is under the threshold
                if e_value < e_value_threshold:
                    if is_valid_description(description):
                        if e_value < min_e_value:
                            min_e_value = e_value
                            corresponding_description = description
                            corresponding_colid = col
                    elif e_value < fallback_min_e_value:  # Track fallback option
                        fallback_min_e_value = e_value
                        fallback_description = description
                        fallback_colid = col

        # Use the best valid description found, or fallback if none were valid
        if corresponding_description is None and fallback_description is not None:
            corresponding_description = fallback_description
            corresponding_colid = fallback_colid
            min_e_value = fallback_min_e_value
        elif corresponding_description is None and fallback_description is None:
            corresponding_description = "No Matches to Sequence Found"
            corresponding_colid = "NA"
            min_e_value = "NA"

        # Store the result for the current row
        results.append({
            'Row Index': index,
            'First Column Value': row[df.columns[0]],
            'Minimum E-value': min_e_value if min_e_value != float('inf') else "NA",
            'Column ID': corresponding_colid,
            'Description': corresponding_description
        })
        #count number of invalid descriptions
        new_invalid_count = sum(not is_valid_description(result['Description']) for result in results)

    # Calculate the number of entries replaced and improvement
    if replace == 'force':
        print('Force mode improvement stats: ', flush=True)
        hyp_replaced = total_hyp_entry - new_invalid_count
        improvement = 0
        if total_hyp_entry > 0 and hyp_replaced >= 0:
            improvement = round((hyp_replaced / total_hyp_entry) * 100, 1)
        if hyp_replaced<0:
            print('Note: -replace force can have negative values for improvement. This generally is in situations where the original annotation is also a placeholder annoation that is not included in placeholders.txt but the best forced  replacement value WAS included in placeholders.txt. Using -replace gentle will not attempt to replace original  annotations that do not appear in the placeholders.txt list. Alternatively, you can update the list to include  new placeholder annotations that may appear in your original annotation.\n', flush=True)
        print('Total "' + '", "'.join(term for term in hyp_syn_list) + '" :' + str(total_hyp_entry), flush=True)
        print('Total replaced: ' + str(hyp_replaced), flush=True)
        print('Percent Improvement (%): ' + str(improvement), flush=True)
        print("\n###################################\n", flush=True)
   


    print('Generating list of best hits', flush=True)
    if '-best_hits' in inputs:
        with open(besthit_output, 'w') as bh_out:
            if results:
                for result in results:
                    bh_out.write(f"Row Index: {result['Row Index']}, First Column Value: {result['First Column Value']}, Minimum E-value: {result['Minimum E-value']}, Column ID: {result['Column ID'].strip(' E-value')}, Description: {result['Description']}" + '\n')
            else:
                bh_out('No matching rows found with E-value < ' + str(e_value_threshold)+'\n')

    print('Generating new annotations', flush=True)
    if replace == 'force':
        fasta_file_path = cds_path  

        base_name = os.path.splitext(os.path.basename(fasta_file_path))[0]  
        new_fasta_file = os.path.join(out_dir_path, f"{base_name}.updated.faa") 

        
        with open(fasta_file_path, 'r') as fasta_file:
            fasta_lines = fasta_file.readlines()

        updated_fasta_lines = []
        for line in fasta_lines:
            if line.startswith('>'):  
                for result in results:
                    first_column_value = result['First Column Value']
                    if first_column_value in line:  
                        line = '>' + str(first_column_value) + ' ' + str(result['Description']) + '\n'
                        break  
            updated_fasta_lines.append(line)

        with open(new_fasta_file, 'w') as new_file:
            new_file.writelines(updated_fasta_lines)

        print(f"Updated FASTA file written to: {new_fasta_file}", flush=True)
        print("\n###################################\n", flush=True)
##################################
#clean databases
##################################
if clean == 'TRUE':
    for name in hmmer_todo:
        shutil.rmtree(out_dir_path + '/' + name)
    for name in blast_todo:
        shutil.rmtree(out_dir_path + '/' + name)
##################################

end = time.time()
print('Time elapsed: ' + str(round((end - start) / 60, 2)) + ' minutes, ' + str(round(((end - start) / 3600), 3)) + ' hours', flush=True)

