#from __future__ import print_function
from datetime import datetime
#import yaml  #CiteSoftLocal will assume these are not available.
#import semantic_version #CiteSoftLocal will assume these are not available.
import re
import sys
import os

def eprint(*args, **kwargs):#Print to stderr
    print(*args, file=sys.stderr, **kwargs)

######GLOBAL VARIABLES SECTION######
citations_dict = {}
checkpoint_log_filename = "CiteSoftwareCheckpointsLog.txt"
consolidated_log_filename = "CiteSoftwareConsolidatedLog.txt"
uniqueIDs_log_filename  = "CiteSoftwareUniqueIDs.txt"
consolidated_CFF_filename = "CITATIONS.cff"
validate_on_fly = True#Flag.  If set to true, argument names will be checked in real time, and invalid argument names will result in a printed warning to the user
valid_optional_fields = ["version", "cite", "author", "doi", "url", "encoding", "misc"]
valid_required_fields = ['timestamp', 'unique_id', 'software_name']
######END OF GLOBAL VARIABLES SECTION######

#TODO: Make a global (in the module) for the absolute path to the log.  That way other modules can change the file_path through CiteSoft.set_log_abs_path(...) or by CiteSoft.log_abs_path = ...
#cwd = os...path....
#log_abs_path = cwd
# def set_log_abs_path(abs_path)
    # global log_abs_path
    # log_abs_path = abs_path

########## START SECTION OF CODE WHICH HAS FUNCTIONS THAT ARE TO BE USED AS DECORATORS ################
#The below functions are intended to be used as decorators.
#They are similar to the example "decorator_maker_with_arguments" at https://www.datacamp.com/community/tutorials/decorators-python
#To find the example, search for "decorator_maker_with_arguments" at the above link.
#function "inner" below is named 'decorator' in the above link and 'wrapper' below is named 'wrapper' in the above link.


#The function_call_cite function is intended to be used for creating a citation any time the wrapped function is called.
#This is expected to be the most common usage of CiteSoft.
#write_immediately means that the citations would be written immediately when encountered, not just stored in memory.
def function_call_cite(unique_id, software_name, write_immediately=True, **add_args):
    #the unique_id and the software_name are the only truly required args.
    #Optional args are: ["version", "cite", "author", "doi", "url", "encoding", "misc"]
    #Every arg must be a string.
    def inner(func):
        def wrapper(*args, **kwargs):
            add_citation(unique_id, software_name, write_immediately, **add_args)
            result = func(*args, **kwargs)
            return result
        return wrapper
    return inner

#function_call_cite used to be called module_call_cite, so the below line is for backwards compatibility.
module_call_cite = function_call_cite

#The after_call_compile_checkpoints_log function is intended to be used as a decorator and only adds citations *after* the function call is finished.
def after_call_compile_checkpoints_log(file_path="", empty_checkpoints=True):
    def inner(func):
        def wrapper(*args, **kwargs):
            result = func(*args, **kwargs)
            compile_checkpoints_log(file_path=file_path, empty_checkpoints=empty_checkpoints)
            return result
        return wrapper
    return inner

#The after_call_compile_consolidated_log function is intended to be used as a decorator and makes the consolided log *after* the function call is finished.
def after_call_compile_consolidated_log(file_path="", compile_checkpoints=False):
    def inner(func):
        def wrapper(*args, **kwargs):
            result = func(*args, **kwargs)
            compile_consolidated_log(file_path=file_path, compile_checkpoints=compile_checkpoints)
            return result
        return wrapper
    return inner

########## END SECTION OF CODE WHICH HAS FUNCTIONS THAT ARE TO BE USED AS DECORATORS ################

#The import_cite function is intended to be used at the top of a sofware module. This would be used when a module should be cited each time it is imported.
#write_immediately means that the citations would be written immediately when encountered, not just stored in memory.
def import_cite(unique_id, software_name, write_immediately=True, **kwargs):
    add_citation(unique_id, software_name, write_immediately, **kwargs)

#The add_citation function is the method which actually adds a citation.
#It's not intended to be called directly by the dev-user, but it could be.
def add_citation(unique_id, software_name, write_immediately=True, **kwargs):
    new_entry = {'unique_id' : unique_id, 'software_name' : software_name, 'timestamp' : get_timestamp()}
    for key in kwargs:
        if validate_on_fly:
            if not key in valid_optional_fields:
                eprint("Warning, " + key + " is not an officially supported argument name.  Use of alternative argument names is strongly discouraged.")
        if type(kwargs[key]) is not list:#Make sure single optional args are wrapped in a list
            kwargs[key] = [kwargs[key]]
        new_entry[key] = kwargs[key]
    if unique_id in citations_dict:#Check for duplicate entries(e.g. from calling the same function twice)
        citations_dict[unique_id] = compare_same_id(citations_dict[unique_id], new_entry)
    else:
        citations_dict[unique_id] = new_entry
    if write_immediately == True:
        compile_checkpoints_log()

#This function creates cff files for each entry based. Th CFF files name is the unique_id converted to a valid file name.
def create_cff(entry, file_path=""):
    if "CITATIONS" not in os.listdir():
        os.mkdir("./CITATIONS")
    import re
    unique_id_string = entry['unique_id']
    valid_file_name_string = re.sub('[^\w_.)( -]', '_', unique_id_string)#remove characters from unique_id_string that are disallowed in filenames strings. Replace with "_". #TODO: We should using an encoding rather than this simple replace, that way the unique_id conversion will be 1:1 and reversible. The reason I did not use an encoding (yet) is I wanted to make sure that the encoding used is easily replicated in all programming languages. It should not be python specific -- but maybe it would be okay if it is?
    cff_filename = valid_file_name_string + ".cff"
    if os.path.exists("./" + file_path + "/CITATIONS/IndividualCff/"):
        pass
    else:
       os.mkdir("./" + file_path + "/CITATIONS/IndividualCff/")
    with open("./" + file_path + "/CITATIONS/IndividualCff/" + cff_filename, 'w') as file:
        write_dict_to_cff(file, entry)
    #now write to the consolided CFF file (since we are appending, it will be created if it does not exist):
    with open("./" + file_path + "/CITATIONS/" + consolidated_CFF_filename, 'a') as file:
        file.write('---\n')
        write_dict_to_cff(file, entry)

#As of August 4th 2021, we will just make minimal CFF files with the minimal required fields:
# https://github.com/citation-file-format/citation-file-format/blob/main/examples/1.2.0/pass/minimal/CITATION.cff
def write_dict_to_cff(file, citation_dict):
    file.write('cff-version: 1.2.0\n')
    file.write('message:' + citation_dict['timestamp'] + '\n')
    file.write('title: ' + citation_dict['software_name'] + '\n')
    file.write('identifiers: \n')
    file.write(' - type: unique_id' + '\n')
    file.write('   value: '+ citation_dict['unique_id'] + '\n')
    if 'author' in citation_dict: #TODO: change this to have family-names and given-names. Also add orcid.
        file.write('authors: \n') 
        for authorName in citation_dict['author']:
            file.write(' - name: ' + authorName + '\n')
    other_cff_valid_fields = ["version", "doi", "url", "license"]
    for field in other_cff_valid_fields:
        if field in citation_dict:
            file.write(field + ": " + str(citation_dict[field][0]) + '\n') #Consider changing: currently CiteSoft makes all optional fields into a list, including the version number. So we are taking the first item in the list. 

def update_unique_IDs_file(file_path=""):
    if "CITATIONS" not in os.listdir():
        os.mkdir("./CITATIONS")
    written_unique_IDs = []
    if uniqueIDs_log_filename in os.listdir("./CITATIONS"): #check if the file exists already.
        #if it exists, grab the items from it and put them into written_unique_IDs.
        with open("./" + file_path + "/CITATIONS/" + uniqueIDs_log_filename, 'r') as file:
            written_unique_IDs = file.readlines()            
    #To consider: Do we need to specify what kind of linebreak characters are used? No, just use "\n" in the specifications.
    for unique_id_index, unique_id in enumerate(written_unique_IDs):
        written_unique_IDs[unique_id_index] = unique_id[:-1] #remove the \n at the end of each string.
    
    #remove the existing unique_IDs from the citations_dict.
    current_keys = list(citations_dict.keys())
    for dict_key in current_keys:
        current_unique_id = citations_dict[dict_key]['unique_id']
        if current_unique_id in written_unique_IDs:
            del citations_dict[current_unique_id]
    
    #write the remaining unique_id values.
    with open("./" + file_path + "/CITATIONS/" + uniqueIDs_log_filename, 'a') as file:
        for dict_key in citations_dict:
            current_unique_id = citations_dict[dict_key]['unique_id']
            if current_unique_id not in written_unique_IDs:
                file.write(current_unique_id + "\n")

#Normally, checkpoints are stored in a dictionary until they are exported.  
#The exporting happens either when requested to from add_citation or from compile_consolidated_log.
def compile_checkpoints_log(file_path="", empty_checkpoints=True):
    if "CITATIONS" not in os.listdir():
        os.mkdir("./CITATIONS")
    update_unique_IDs_file(file_path=file_path) #write any unique ideas
    with open("./" + file_path + "/CITATIONS/" + checkpoint_log_filename, 'a') as file:
        write_dict_to_output(file, citations_dict)
    for dict_key in citations_dict:
        create_cff(citations_dict[dict_key])
    if empty_checkpoints==True:
        citations_dict.clear()

def compile_consolidated_log(file_path="", compile_checkpoints=False):
    #CiteSoftLocal will use compile checkpoints since it cannot consolidate.
    compile_checkpoints_log()
    print("Warning: CiteSoftLocal cannot make a consolidated log. Citations have been exported to CiteSoftwareCheckpointsLog.txt")
    # ##BELOW IS THE ORIGINAL FUNCTION CODE3##
    # #compile_checkpoints =True must be used sparingly. Otherwise it can slow down performance when the checkpoints file gets large.
    # import copy
    # consolidated_dict = copy.deepcopy(citations_dict) #Make a copy of citations_dict before emptying it.
    # if compile_checkpoints == True:
        # compile_checkpoints_log()
    # #Grab the citations from the checkpoint log.        
    # if checkpoint_log_filename in os.listdir("./CITATIONS"): #check if the file exists already.
        # checkpoint_log_exists = True
    # else:
        # checkpoint_log_exists = False
    # if checkpoint_log_exists == True: #can only read file if it exists.
        # with open("./CITATIONS/" + checkpoint_log_filename, 'r') as file:
            # yaml_file_contents = yaml.safe_load_all(file)
            # for yaml_document in yaml_file_contents:
                # if yaml_document != None: #This is for 'blank' documents of "---" with nothing after that symbol.
                    # for citation_entry in yaml_document:
                        # id = citation_entry["unique_id"]
                        # if id in consolidated_dict:
                            # consolidated_dict[id] = compare_same_id(consolidated_dict[id], citation_entry)
                        # else:
                            # consolidated_dict[id] = citation_entry        
    # #Grab the citations from the consolidated log.    
    # if consolidated_log_filename in os.listdir("./CITATIONS"): #check if the file exists already.
        # consolidated_log_exists = True
    # else:
        # consolidated_log_exists = False
    # if consolidated_log_exists == True: #can only read file if it exists.
        # with open("./" + file_path + "/CITATIONS/" + consolidated_log_filename, "r") as file:
            # yaml_file_contents = yaml.safe_load_all(file)
            # for yaml_document in yaml_file_contents:
                # if yaml_document != None: #This is for 'blank' documents of "---" with nothing after that symbol.
                    # for citation_entry in yaml_document:
                        # id = citation_entry["unique_id"]
                        # if id in consolidated_dict:
                            # consolidated_dict[id] = compare_same_id(consolidated_dict[id], citation_entry)
                        # else:
                            # consolidated_dict[id] = citation_entry
    # with open("./CITATIONS/" + consolidated_log_filename, 'w') as file:
        # file.write('#Warning: CiteSoftwareConsolidatedLog.txt may not include all softwares used. It is the end-user’s responsibility to verify that the no software citations are missing relative to those recorded in the complete logfile, CiteSoftwareCheckpointsLog.txt . This verification is important to do when using two or more codes together for the first time.\n')
        # write_dict_to_output(file, consolidated_dict)

#Takes a dictionary, converts it to CiteSoft-compatible YAML, and writes it to file
def write_dict_to_output(file, dictionary):
    if len(dictionary) > 0:
        file.write('---\n')
        for key,dict in dictionary.items():
            file.write('-\n')
            for s in valid_required_fields:
                file.write('    ' + s + ': >-\n')
                file.write('    '*2 + dict[s] + '\n')
            for subkey in dict:
                if subkey not in valid_required_fields:
                    file.write('    ' + subkey + ':\n')
                    if type(dict[subkey]) is list:
                        for i in dict[subkey]:
                            file.write('    '*2 + '- >-\n')
                            file.write('    '*3 + i + '\n')
                    else:
                        file.write('    '*2 + '- >-\n')
                        file.write('    '*3 + dict[subkey] + '\n')

#Helper Functions

#Returns a string of the current time in the ISO 8601 format (YYYY-MM-DDThh:mm:ss).
def get_timestamp():
    now = datetime.now()
    timestamp = now.strftime("%Y-%m-%dT%H:%M:%S")
    return timestamp

#Compares two entries
#Returns : The entry which should be kept
def compare_same_id(old_entry, new_entry):
    return new_entry #CiteSoftLocal will not do comparisons. It will just return the new_entry.
    # ##BELOW IS THE ORIGINAL FUNCTION CODE##
    # old_has_version = "version" in old_entry
    # new_has_version = "version" in new_entry
    # if old_has_version and new_has_version:#If both entries have a version, compare them return and the return the greater(newer) version
        # old_ver_str = str(old_entry["version"][0])
        # new_ver_str = str(new_entry["version"][0])
        # #Initialize variables, assume strings are valid unless parsing fails
        # old_ver_semver_valid = True
        # new_ver_semver_valid = True
        # decimal_regex_str = "^[0-9]+\.[0-9]+$"#Match string with decimal point enclosed by at least one number on either side
        # if re.match(decimal_regex_str, old_ver_str):
            # old_ver_str += '.0'#To ensure semantic version parser handles a decimal value correctly
        # if re.match(decimal_regex_str, new_ver_str):
            # new_ver_str += '.0'#To ensure semantic version parser handles a decimal value correctly
        # try:
            # old_sv = semantic_version.Version(old_ver_str)
        # except ValueError:
            # old_ver_semver_valid = False
        # try:
            # new_sv = semantic_version.Version(new_ver_str)
        # except:
            # new_ver_semver_valid = False
        # if old_ver_semver_valid and new_ver_semver_valid:#If both entries have a valid SemVer version, keep the older one only if it's greater. Else, keep the newer one.
            # if old_sv > new_sv:
                # return old_entry
            # else:
                # return new_entry
        # elif old_ver_semver_valid:#If only the old entry has a valid SemVer version, keep it
            # return old_entry
        # elif new_ver_semver_valid:#If only the new entry has a valid SemVer version, keep it
            # return new_entry
        # else:
            # #Version comparison failed, use alphanumeric comparison
            # if old_ver_str > new_ver_str:
                # return old_entry
            # else:
                # return new_entry
    # elif old_has_version and not new_has_version:#If old entry has a version and the new entry doesn't, the entry with a version takes precedence
        # return old_entry
    # elif not old_has_version and new_has_version:#Likewise, if new entry has a version and the old entry doesn't, the entry with a version takes precedence
        # return new_entry
    # else:#If neither entry has a version, save the new entry
        # return new_entry
