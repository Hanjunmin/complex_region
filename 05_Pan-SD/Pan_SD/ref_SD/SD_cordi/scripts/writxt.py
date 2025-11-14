
import os
import glob
import pickle
import threading
import multiprocessing

def write_pairset(file_name, pairset):
    with open(file_name + '.a.txt', 'w') as file:
        for key, value in pairset.items():
            for val in value:
                file.write(f"{key}\t{val}\n")

def write_jacset(file_name, jacset):
    with open(file_name + '.b.txt', 'w') as file:
        for key, value in jacset.items():
            for val in value:
                file.write(f"{key}\t{val}\n")

def write_otherid(file_name, otherid):
    with open(file_name + '.c.txt', 'w') as file:
        for value in otherid:
            file.write(f"{value}\n")

def process_file(file):
    file_name = os.path.splitext(os.path.basename(file))[0]
    with open(file, 'rb') as f:
        loaded_vars = pickle.load(f)
    print(file_name)
    pairset = loaded_vars['var2']
    jacset = loaded_vars['var3']
    otherid = loaded_vars['var1']
    loaded_vars.clear()
    
    # Start threads for writing files
    pairset_thread = threading.Thread(target=write_pairset, args=(file_name, pairset))
    jacset_thread = threading.Thread(target=write_jacset, args=(file_name, jacset))
    otherid_thread = threading.Thread(target=write_otherid, args=(file_name, otherid))
    pairset_thread.start()
    jacset_thread.start()
    otherid_thread.start()
    pairset_thread.join()
    jacset_thread.join()
    otherid_thread.join()
    
    pairset.clear()
    jacset.clear()
    otherid.clear()


current_directory = os.getcwd()
files = glob.glob(os.path.join(current_directory, '*end.pkl'))
with multiprocessing.Pool(processes=30) as pool:
    pool.map(process_file, files)

