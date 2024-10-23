print('Loading libraries...')


'''STEPS FOR RUNNING FILE:

            - Create directory for system of interest
            - Run getxyz.py to make folders required for runs
            - ensure that in the directory above the working directory energy.inp and minmise.inp are included
            - Change the job script template for each system, and change BASIS SET and GTH_POTENTIALS file location for each system
            - Run script and follow instructions'''

#differing strain jiggle and vol options can be used, alter here

vol_list = [-10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110,
                        120, 130, 140, 150, 160, 170, 180, 190, 200] #different volume supercells
strain_list = (-5, -2, 2, 5)  #different strain percentages for each volume
n_jiggle = 5 #number of different perturbations performed for each strained structure
n_strain = 5


import numpy as np
from pymatgen.core import Structure
from pymatgen.transformations.standard_transformations import PerturbStructureTransformation, DeformStructureTransformation
from pymatgen.io.xyz import XYZ
import os
import shutil
import sys
import time
import subprocess
import ast 

print('libraries loaded!')


'''TO DO:
        - decide on number of atoms for CELL optimisation
        - change number of requested nodes in job submission
    '''

def vol_scale(file_name, scale_factor, out_file_name):

    coords, atoms, a, b, c = read_xyz(file_name)
    lattice = lattice_vectors(file_name)

    structure = Structure(lattice, atoms, coords, coords_are_cartesian = True)

    volume_scale_factor = float(scale_factor) ** (1 / 3)

    new_lattice = structure.lattice.matrix * volume_scale_factor
    new_coords = coords * volume_scale_factor

    scaled_structure = Structure(new_lattice, atoms, new_coords, coords_are_cartesian=True)

    a, b, c = scaled_structure.lattice.abc
    alpha, beta, gamma = scaled_structure.lattice.angles

    lattice2 = scaled_structure.lattice.matrix
    vector_string = f'[[{lattice2[0][0]},{lattice2[0][1]},{lattice2[0][2]}],[{lattice2[1][0]},{lattice2[1][1]},{lattice2[1][2]}],[{lattice2[2][0]},{lattice2[2][1]},{lattice2[2][2]}]]'
    
    xyz = XYZ(scaled_structure)
    xyz_string = xyz.__str__()

    #add lattice parameters and angles as comment
    comment = f"Lattice: a={a:.6f} b={b:.6f} c={c:.6f}; Angles: alpha={alpha:.6f} beta={beta:.6f} gamma={gamma:.6f}; Vectors: {vector_string}"
    xyz_lines = xyz_string.split('\n')
    xyz_lines.pop(1) #Remove original comment line
    xyz_lines.insert(1, comment)
    modified_xyz_string = '\n'.join(xyz_lines)
    #write to XYZ file
    with open(out_file_name, "w") as f:
        f.write(modified_xyz_string)

   

    return out_file_name

def replace_line(filename, string_to_replace, new_string):
        with open(filename, 'r') as f:
            lines = f.readlines()
        with open(filename, 'w') as f:
            for line in lines:
                if string_to_replace in line:
                    line = line.replace(string_to_replace, new_string)
                    f.write(line)
                else:
                    f.write(line)


def read_xyz(xyz_file):
        coords = []
        atoms = []
        with open(xyz_file, 'r') as inp:
            lines = inp.readlines()
            comment = lines[1].split(';')[0].split(':')[1].split()
            a = float(comment[0].split('=')[1])
            b = float(comment[1].split('=')[1])
            c = float(comment[2].split('=')[1])

            for line in lines[2:]:
                if line.strip():  # Ensure non-empty line
                    atom, x, y, z = line.split()[:4]
                    atoms.append(atom)
                    coords.append([float(x), float(y), float(z)])

        coords = np.array(coords)
        atoms = np.array(atoms)
        return coords, atoms, a, b, c

def lattice_vectors (file):
    with open(file, 'r') as f:
          vectors = f.readlines()[1].split(';')[-1].split(':')[1].strip()
          lattice_vectors = ast.literal_eval(vectors)
    return lattice_vectors



def strain(input_file, out_file_name, strain_range = (-0.05, 0.05)):
    coords, atoms, a, b, c = read_xyz(input_file)
    lattice = lattice_vectors(input_file)

    #create structure object
    structure = Structure(lattice, atoms, coords, coords_are_cartesian = True)

    displacements = [np.random.uniform(*strain_range) for _ in range(3)]
    zero_index = np.random.randint(0, 3)  # One of the displacements is 0
    displacements[zero_index] = 0




    deformation_matrix = [[1 + displacements[0], 0, 0], [0, 1 + displacements[1], 0], [0, 0, 1 + displacements[2]]]
    strain = DeformStructureTransformation(deformation_matrix)
    deformed_xyz = strain.apply_transformation(structure)
    lattice2 = deformed_xyz.lattice.matrix
    vector_string = f'[[{lattice2[0][0]},{lattice2[0][1]},{lattice2[0][2]}],[{lattice2[1][0]},{lattice2[1][1]},{lattice2[1][2]}],[{lattice2[2][0]},{lattice2[2][1]},{lattice2[2][2]}]]'


    a, b, c = deformed_xyz.lattice.abc
    alpha, beta, gamma = deformed_xyz.lattice.angles

    xyz = XYZ(deformed_xyz)
    xyz_string = xyz.__str__()


    #add lattice parameters and angles as comment
    comment = f"Lattice: a={a:.6f} b={b:.6f} c={c:.6f}; Angles: alpha={alpha:.6f} beta={beta:.6f} gamma={gamma:.6f}; Vectors: {vector_string}"
    xyz_lines = xyz_string.split('\n')
    xyz_lines.pop(1) #Remove original comment line
    xyz_lines.insert(1, comment)
    modified_xyz_string = '\n'.join(xyz_lines)
    #write to XYZ file
    with open(out_file_name, "w") as f:
        f.write(modified_xyz_string)



def jiggle (input_file, out_file_name, max_min = (0.5, 0.05)):

    coords, atoms, a, b, c = read_xyz(input_file)
    lattice = lattice_vectors(input_file)

    
    #create structure object
    structure = Structure(lattice, atoms, coords, coords_are_cartesian = True)

    #Petrubation matrix (between 0.05 and 0.5 angstrom deformation)
    perturbation = PerturbStructureTransformation(*max_min)


    perturbed_structure = perturbation.apply_transformation(structure)

    a, b, c = perturbed_structure.lattice.abc
    alpha, beta, gamma = perturbed_structure.lattice.angles
    lattice2 = perturbed_structure.lattice.matrix
    vector_string = f'[[{lattice2[0][0]},{lattice2[0][1]},{lattice2[0][2]}],[{lattice2[1][0]},{lattice2[1][1]},{lattice2[1][2]}],[{lattice2[2][0]},{lattice2[2][1]},{lattice2[2][2]}]]'

    xyz = XYZ(perturbed_structure)
    xyz_string = xyz.__str__()

    #add lattice parameters and angles as comment
    comment = f"Lattice: a={a:.6f} b={b:.6f} c={c:.6f}; Angles: alpha={alpha:.6f} beta={beta:.6f} gamma={gamma:.6f}; Vectors: {vector_string}"
    xyz_lines = xyz_string.split('\n')
    xyz_lines.pop(1) #Remove original comment line
    xyz_lines.insert(1, comment)
    modified_xyz_string = '\n'.join(xyz_lines)
    #write to XYZ file

    with open(out_file_name, "w") as f:
        f.write(modified_xyz_string)



def supercell (input_file, x, y, z):
    coords, atoms, a, b, c = read_xyz(input_file)
    lattice = lattice_vectors(input_file)
   
    #create structure object
    structure = Structure(lattice, atoms, coords, coords_are_cartesian = True)
    supercell = structure.make_supercell([[x, 0, 0], [0, y, 0], [0, 0, z]])

    a, b, c = supercell.lattice.abc
    alpha, beta, gamma = supercell.lattice.angles
    lattice2 = supercell.lattice.matrix
    vector_string = f'[[{lattice2[0][0]},{lattice2[0][1]},{lattice2[0][2]}],[{lattice2[1][0]},{lattice2[1][1]},{lattice2[1][2]}],[{lattice2[2][0]},{lattice2[2][1]},{lattice2[2][2]}]]'

    xyz = XYZ(supercell)
    xyz_string = xyz.__str__()
    num_atoms = int(xyz_string.splitlines()[0])


    #add lattice parameters and angles as comment
    comment = f"Lattice: a={a:.6f} b={b:.6f} c={c:.6f}; Angles: alpha={alpha:.6f} beta={beta:.6f} gamma={gamma:.6f}; Vectors: {vector_string}"
    xyz_lines = xyz_string.split('\n')
    xyz_lines.pop(1) #Remove original comment line
    xyz_lines.insert(1, comment)
    modified_xyz_string = '\n'.join(xyz_lines)
    return modified_xyz_string, num_atoms


def extract_comment_change_inp (input_file, xyz_file, symmetry = None):
    with open(xyz_file, 'r') as f:
        comment = f.readlines()[1]
        x, y, z = comment.split(';')[0].split(':')[1].split()[0:3]
        ABC_string = f"{x.split('=')[1]} {y.split('=')[1]} {z.split('=')[1].replace(';', '')}"
        a, b, c = comment.split(';')[1].split(':')[1].split()[0:3]
        abc_string = f"{a.split('=')[1]} {b.split('=')[1]} {c.split('=')[1].replace(';', '')}"

    with open(input_file, 'r') as f:
        lines = f.readlines()
        for index, line in enumerate(lines):
            if '$ABC' in line:
                lines[index] = line.replace('$ABC', ABC_string)
            if '$abc' in line:
                lines[index] = line.replace('$abc', abc_string)

    with open(input_file, 'w') as f:
        f.writelines(lines)



def check_completion(start_index):
    current_vol = None
    current_strain = None
    current_jiggle = None
    with open("progress.txt", "r") as f:
        lines = f.readlines()
        for index, line in enumerate(lines):

            line = line.strip()

            if not line or line.startswith('+'):
                continue

            columns = line.split('|')
            columns = [col.strip() for col in columns]

            if columns[1].startswith("vol"):
                current_vol = columns[1]
                current_strain = None
                current_jiggle = None

            elif columns[2].startswith("strain"):
                current_strain = columns[2]
                current_jiggle = None

            elif columns[3].startswith("jiggle"):
                current_jiggle = columns[3]
            # Fetch path
            if current_jiggle is not None:
                path = f'{current_vol}/{current_strain}/{current_jiggle}/'
            elif current_strain is not None:
                path = f'{current_vol}/{current_strain}/'
            elif current_vol is not None:
                path = f'{current_vol}/'

            if index > start_index:
                simulation_complete_check = False
                # If simulation not yet run
                if columns[4] == 'n':
                    lines[index] = lines[index].replace(' n ', ' r ')   #when running sim, complete column is changed to r, for running
                    break

                # If simulation has been run
                if columns[4] == 'r':

                    #navigate to directory where simulation has been run
                    path2 = os.path.join(path, 'result.out')
                    if os.path.exists(path2): #if simulation has started
                        with open (path2, 'r') as f2:
                            lines2 = f2.readlines()
                        #Determine whether simulation has been completed
                        for line2 in lines2:
                            if 'The number of warnings for this run is' in line2:
                                simulation_complete_check = True
                        if simulation_complete_check == False: #simulation has not been complete but has started.

                            # Check job is running
                            result = subprocess.run(['squeue', '--me', '--format="%j"'],
                                                capture_output = True, text = True)
                            if result.returncode == 0:
                                job_names = [line for line in result.stdout.splitlines() if line]
                                for job in job_names:
                                    job = job.strip('"').strip()
                                    path = path.strip()
                                    if str(job) == str(path):
                                        #if simulation is in queue we skip running it
                                        print(f'Job {path} already queued, skipping...')
                                        simulation_complete_check = True
                                        break
                                #Only executes if simulation is not in queue
                                # i.e. if simulation has started but not finished and is not in queue, we restart
                                print(f'ERROR: {path} has started but is not currently running. Restarting simulation (CHECK THERE IS ENOUGH TIME TO COMPLETE SIMULATION).')
                                break

                            else:
                                print('error in checking queue status')


                        else: #if simluation complete check returns true
                            lines[index] = lines[index].replace(' r ', ' y ')




                    else: #simulation has not started, or is still in queue
                        #check if simulation is in queue
                        result = subprocess.run(['squeue', '--me', '--format="%j"'],
                                                capture_output = True, text = True)
                        if result.returncode == 0:
                            job_names = [line for line in result.stdout.splitlines() if line]
                            for job in job_names:
                                if job == path:
                                    #if simulation is in queue and not started yet we skip running it
                                    print('Job already queued, skipping...')
                                    simulation_complete_check = True
                                    break
                            if simulation_complete_check == False:
                                #if simulation is not in queue and has not started yet but should have run
                                print(f'ERROR: simulation {path} should have been run but has not been. Restarting...')
                                break


                        else:
                            print('error in checking queue status')


    with open ('progress.txt', 'w') as o:
        o.writelines(lines)
    print(f'Submitting job {path}')


    return current_vol, current_strain, current_jiggle, simulation_complete_check, index



def create_job_script (name, input_path, output_path, tag, time, nodes):
    job_script = f"""#!/bin/bash

#SBATCH --job-name={tag}
#SBATCH --nodes={nodes}
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=8
#SBATCH --time={time}

#SBATCH --account=e05-power-mur
#SBATCH --partition=standard
#SBATCH --qos=standard


source /work/e05/e05/ash141/myenv/bin/activate
module load PrgEnv-gnu
module load cray-fftw
module load cray-python
module load mkl

#Ensure OMP_NUM_THREADS is consistent with cpus-per-task above
export OMP_NUM_THREADS=8
export OMP_PLACES=cores
c
export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK

srun --hint=nomultithread --distribution=block:block /work/e05/e05/ash141/cp2k/cp2k/cp2k-2024.1/exe/ARCHER2/cp2k.psmp -i {input_path} -o {output_path}
"""

    with open (f'{name}.slurm', 'w') as inp:   #Create submission script
        inp.write(job_script)


def run_simulation (start_index):
    current_vol, current_strain, current_jiggle, simulation_complete_check, index = check_completion(start_index)  #Already running sims have r in last column, check for next sim with n

    if simulation_complete_check == False:
    #Determine path to inp file that is scheduled to run
        if current_jiggle is not None:
            path = f'{current_vol}/{current_strain}/{current_jiggle}/'
        elif current_strain is not None:
            path = f'{current_vol}/{current_strain}/'
        elif current_vol is not None:
            path = f'{current_vol}/'

        os.chdir(path)

        #NEED TO CHANGE SPECIFICS DEPENDING ON SYSTEM, CHECK JOB SCRIPT STRUCTURE SO IT APPLIES TO HPC
        create_job_script(name = 'subjob', input_path= f'energy.inp', output_path= f'result.out', tag = f'{path}', time = '01:00:00', nodes = '4')

        os.system('sbatch subjob.slurm')   #Submit new job to queue
        os.chdir(parent)

    return current_vol, current_strain, current_jiggle, simulation_complete_check, index








if __name__ == "__main__":

    parent = os.getcwd()




#################################### OPTIMISE SIZE OF SUPERCELL ##############################################################################
    path = os.path.join(parent, 'Cells')

    if not os.path.exists(path):


        coords, atoms, a, b, c = read_xyz('inp.xyz')
        #Find cell multipliers (supercell size is proportional to unit cell size for equal sampling in k space)
        largest_param = np.max([float(a), float(b), float(c)])
        a_multiplier = round(largest_param / a)
        b_multiplier = round(largest_param / b)
        c_multiplier = round(largest_param / c)

        shutil.copy('../energy.inp', '.')

        os.mkdir('Cells')
        os.chdir('Cells')
        os.mkdir('1')
        os.chdir('1')
        shutil.copy('../../inp.xyz', '.')
        shutil.copy('../../energy.inp', '.')
        extract_comment_change_inp(input_file= 'energy.inp', xyz_file= 'inp.xyz')

        with open('energy.inp', 'r') as f:
            lines = f.readlines()
            for index, line in enumerate(lines):
                if 'supercell_min.xyz' in line:
                    lines[index] = line.replace('supercell_min.xyz', 'inp.xyz')

        with open('energy.inp', 'w') as f:
            f.writelines(lines)

        create_job_script(name= 'job', input_path= 'energy.inp', output_path= 'result.out', tag = 'CELL', time = '01:00:00', nodes = '4')

        os.system('sbatch job.slurm')
        print('Submited 1x1x1 cell!')

        os.chdir('../')
        num_atoms = 0
        N = 1
        while num_atoms < 1000:

            a = (N * a_multiplier)
            b = (N * b_multiplier)
            c = (N * c_multiplier)



            modified_xyz_string, num_atoms = supercell(input_file = '../inp.xyz', x = a, y = b, z = c )

            if num_atoms > 1000:
                break

            os.mkdir(f'{N + 1}')
            os.chdir(f'{N + 1}')
            shutil.copy('../../energy.inp', '.')
            shutil.copy('../../inp.xyz', '.')

            with open('inp.xyz', "w") as f:
                f.write(modified_xyz_string)

            extract_comment_change_inp(input_file= 'energy.inp', xyz_file= 'inp.xyz')

            with open('energy.inp', 'r') as f:
                lines = f.readlines()
                for index, line in enumerate(lines):
                    if 'supercell_min.xyz' in line:
                        lines[index] = line.replace('supercell_min.xyz', 'inp.xyz')

            with open('energy.inp', 'w') as f:
                f.writelines(lines)

            create_job_script(name= 'job', input_path= 'energy.inp', output_path= 'result.out', tag = 'CELL', time = '2:00:00', nodes = '8')

            os.system('sbatch job.slurm')
            print(f'Submited {a}x{b}x{c} cell!')

            os.chdir('../')

            N += 1

        # All simulations are now in queue for CELL optimisation
        print('Once CELL simulations complete, proceed to run minimisation.')
        print('Optimisation running, exiting script!')
        sys.exit()
    else:
        os.chdir('Cells')

        files = os.listdir()
        dirs = [file for file in files if os.path.isdir(file)]

        for dir in range(1, (len(dirs) + 1)):
            path = os.path.join(os.getcwd(), f'{dir}/result.out')
            if os.path.exists(path):
                with open (f'{dir}/result.out', 'r') as f:
                    lines = f.readlines()
                simulation_complete_check = False
                for line in lines:
                    if 'The number of warnings for this run is' in line:
                        simulation_complete_check = True
                if simulation_complete_check == False:
                    print(f'Cell simulations incomplete! (directory {os.getcwd()})\nCheck files and rerun incomplete simulations.')
                    sys.exit()
            else:
                print(f'Cell simulations incomplete! (directory {os.getcwd()})Check files and rerun incomplete simulations.')
                sys.exit()


################################################# PERFORM MINIMSATION ON CONVERGED CELL SIZE #######################################################################


    # Supercell convergence is complete

    files = os.listdir()
    dirs = [file for file in files if os.path.isdir(file)]
    column_widths = [20, 20, 20, 25]


    format_string = "| {{:{w0}}} | {{:{w1}}} | {{:{w2}}} | {{:{w3}}} |\n".format(
        w0=column_widths[0],
        w1=column_widths[1],
        w2=column_widths[2],
        w3=column_widths[3]
    )
    #Write supercell energy data to file
    data_file = os.path.join(os.getcwd(), 'data.txt')
    if not os.path.exists(data_file):
        print('Extracting energies for different supercell sizes and writing to file..')
        with open('data.txt', 'w') as f:
            f.write(format_string.format("Directory", "Energy (a.u.)", "Energy / atom (eV)", "Energy Diff / atom (eV)"))


            f.write(f"{'+' + '-' * (column_widths[0] + 2)}+{'-' * (column_widths[1] + 2)}+{'-' * (column_widths[2] + 2)}+{'-' * (column_widths[3] + 2)}+\n")
            

            # N = 1
            # for dir in range(len(dirs), 0, -1):
            #     dir = str(dir)
                
            #     energy = 0
            #     num_atoms = 0

            #     with open(f'{dir}/result.out', 'r') as r:
            #         lines = r.readlines()
            #         for line in lines:
            #             if 'ENERGY| ' in line:
            #                 energy = float(line.split()[8])
            #             if '- Atoms:' in line:
            #                 num_atoms = int(line.split()[2])

                

            #     energy_per_atom = (energy * 27.2113) / num_atoms if num_atoms > 0 else 0

            #     if N == 1:
            #         E_nmax = energy_per_atom
            #         energy_diff_per_atom = "N/A"

            #     else:
            #         energy_diff_per_atom = E_nmax - energy_per_atom
        

            #     #Convert energy_diff_per_atom to string with formatting, if it's a number
            #     energy_diff_per_atom_str = f"{np.abs(energy_diff_per_atom):.6f}" if energy_diff_per_atom != "N/A" else "N/A"

            #     f.write(format_string.format(
            #         dir,
            #         f"{energy:.6f}",
            #         f"{energy_per_atom:.6f}",
            #         energy_diff_per_atom_str
            #     ))
            #     N += 1
            
            energies = []
            atoms = []
            for dir in range(1, len(dirs) + 1):
                dir = str(dir)
                with open(f'{dir}/result.out', 'r') as r:
                    lines = r.readlines()
                    for line in lines:
                        if 'ENERGY| ' in line:
                            energy = float(line.split()[8])
                            
                        if '- Atoms:' in line:
                            num_atoms = int(line.split()[2])
                    energies.append(energy)
                    atoms.append(num_atoms)

            E_nmax = (float(energies[-1]) * 27.2113) / atoms[-1] #Find energy of largest cell
            for N in range(1, len(energies) + 1):
                energy_per_atom = (energies[N-1] * 27.2113) / atoms[N-1]
                energy_diff_per_atom = E_nmax - energy_per_atom

                energy_diff_per_atom_str = f"{np.abs(energy_diff_per_atom):.6f}" 
                f.write(format_string.format(
                    f"{N}",
                    f"{energy:.6f}",
                    f"{energy_per_atom:.6f}",
                    energy_diff_per_atom_str
                ))


                

    
    min_check = os.path.join(parent, 'MIN.out')
    if not os.path.exists(min_check):
        dir_name = None

        #Read 'data.txt' to find the directory with the smallest energy difference

        with open('data.txt', 'r') as f:
            lines = f.readlines()[2:]
            e_diffs = []
            dir_names = []
            for line in lines:
                e_diff = line.split('|')[4].strip()  
                e_diff_value = float(e_diff) if e_diff != 'N/A' else float('inf')
                e_diffs.append(e_diff_value)
                dir_names.append(line.split('|')[1].strip())

            for i in range(len(e_diffs) - 1):
                if float(e_diffs[i]) < 1E-3 and float(e_diffs[i+1]) < 1E-3: # check that current and subsqeunt energy is within convergence criteria
                    dir_name = dir_names[i]
                    print('Cell within convergence criteria!')
                    break
            if dir_name == None:
                dir_name = dir_names[-1]

        print('Optimal cell size determined')
        #navigate to home
        os.chdir('../')
        #supercell.xyz contains converged cell size
        shutil.copy(f'Cells/{dir_name}/inp.xyz', 'supercell.xyz')
        print('xyz file copied to main directory')


        #supercell dimensions extracted and placed into minimisation file

        shutil.copy(f'../minimise.inp', '.')
        print('Updating minimise.inp with lattice vectors of input cell...')
        extract_comment_change_inp(input_file= 'minimise.inp', xyz_file= 'supercell.xyz')

        print('Minimising cell...')

        #Run Minimisation
        create_job_script(name= 'job', input_path= 'minimise.inp', output_path= 'MIN.out', tag = 'MIN', time = '24:00:00', nodes = '12')
        with open ('inp.xyz', 'r') as f:
            comment_line = f.readlines()[1]
            symmetry = comment_line.split(';')[2].split(':')[1]

        #add in symmetry argument extracted from inp.xyz file
        replace_line('minimise.inp', '$symmetry', symmetry)
        os.system('sbatch job.slurm')   #Submit new job to queue

        print('Submitted minimisation job, resubmit script after minmisation complete to begin deformation simulations.')
        print('Exiting script...')
        sys.exit()

    os.chdir(parent)


    #check if deformation files have already been created
    file_creation_marker = False
    files = os.listdir()
    dirs = [file for file in files if os.path.isdir(file)]
    for directory in dirs:
        base_name = os.path.basename(directory)
        if base_name == 'vol0':
            file_creation_marker = True
            break
    # if not check that minimisation is complete
    if file_creation_marker == False:
        minfile = os.path.join(parent, 'MIN.out')
        min_progress = False
        if os.path.exists(minfile):
            with open('MIN.out', 'r') as f:
                lines = f.readlines()
            for line in lines:
                if 'The number of warnings for this run is' in line:
                    min_progress = True
            if min_progress == False:
                print('Minimisation not yet complete or has an error, rerun.')
                sys.exit()
############################################################ DEFORM MINIMISED CELL, MAKE RELEVENT DIRECTORIES AND WRITE PROGRESS TRACKER #################################


            # for case where files have not yet been created however minimisation is done:

            '''
            1) Previous minimisation spits out positions of cell minimised structure
            2) The min.out file is searched for the optimised cell parameters
            3) These optimised parameters are used to update the lattice vectors of the minimised structure
            4) The new positions and lattice vectors are printed into supercell_min.xyz
            5) The lattice vectors from supercell_min.xyz are parsed foreword into the structure alteration functions
            (vol, strain, jiggle)
            6) The corresponding copy of supercell_min.xyz in the relevant subfolder has its lattice vectors altered
            for the relevant transformation, and the energy.inp file is modified to reflect these new lattice vectors
            7) Generates progress file for checking simulation progress
            '''


            with open('MIN-pos-1.xyz', 'r') as f:
                lines = f.readlines()
                num_atoms = int(lines[0])
                min_structure = lines[-(num_atoms+2):]

            with open('MIN.out', 'r') as f:
                lines = f.readlines()
                for line in lines:
                    if 'CELL| Vector a' in line:
                        a = line.split()[4:7]
                    if 'CELL| Vector b' in line:
                        b = line.split()[4:7]
                    if 'CELL| Vector c' in line:
                        c = line.split()[4:7]
                        print(c)
                lattice_vecs = f'[[{a[0]},{a[1]},{a[2]}],[{b[0]},{b[1]},{b[2]}],[{c[0]},{c[1]},{c[2]}]]'
                print(lattice_vecs)

            with open ('supercell.xyz', 'r') as f:
                comment_line = f.readlines()[1]

            print('Extracted cell parameters from minimisation!')


            with open('supercell_min.xyz', 'w') as f:
                Lattice_comment = comment_line.split(';')[0]
                newline = comment_line.replace(Lattice_comment, f"Lattice: a={float(a[0]):.6f} b={float(b[1]):.6f} c={float(c[2]):.6f}").strip('\n')
                newline = newline + f'; Vectors: {lattice_vecs}\n'
                min_structure[1] = newline
                f.writelines(min_structure)

            print('Written relaxed lattice vectors and positions into supercell_min.xyz!')


            # Calculate lindemann criteria

            coords, atoms, a, b, c,  = read_xyz('inp.xyz')
            lindemann_criteria  = float(min((a, b, c)) * 0.1)


            for N in vol_list:
                dir_name = f'vol{N}'
                if os.path.exists(dir_name) and os.path.isdir(dir_name):
                    print('VOLUME FILE EXISTS ALREADY!')
                    sys.exit()
            
                os.mkdir(f'vol{N}')

            for dir_name in os.listdir('.'):
                if dir_name.startswith('vol') and os.path.isdir(dir_name):
                    vol_path = os.path.join(parent, dir_name)

                    print(f'Deforming for {dir_name}...')
                    if int(dir_name.split('vol')[1]) != 0:
                        scale = float(1 + (float(dir_name.split('vol')[1]))/100)
                    else:
                        scale = 1
                    #take minimised structure (supercell_min.xyz), 
                    vol_scale(file_name= 'supercell_min.xyz', scale_factor= f'{scale:.2f}', out_file_name=f'{vol_path}/supercell_min.xyz') 

                    
                    for s in range(1, (n_strain + 1), 1):
                        xyz = f'strain{s}'
                        strain_path = os.path.join(vol_path, xyz)
                        if os.path.exists(strain_path) and os.path.isdir(strain_path):
                            shutil.rmtree(xyz)
                        os.mkdir(strain_path)
        
                        if s == 1:
                            strain(input_file = f'{vol_path}/supercell_min.xyz', strain_range = (0, 0), out_file_name=f'{strain_path}/supercell_min.xyz')

                        else:
                            strain(input_file = f'{vol_path}/supercell_min.xyz', strain_range = (-0.05, 0.05), out_file_name=f'{strain_path}/supercell_min.xyz')
                        
                       
                        print(f'Deforming for strain {s}...')
                           

                        for N in range(1, (n_jiggle + 1), 1):
                            dir_name = f'jiggle{N}'
                            jiggle_path = os.path.join(strain_path, dir_name)
                            print(f'Performing jiggle {N}...')
                            if os.path.exists(jiggle_path) and os.path.isdir(jiggle_path):
                                shutil.rmtree(dir_name)
                            os.mkdir(jiggle_path)

                            if N ==1:
                                jiggle(f'{strain_path}/supercell_min.xyz', max_min=(0, 0), out_file_name=f'{jiggle_path}/supercell_min.xyz')
                            else:
                                jiggle(f'{strain_path}/supercell_min.xyz', max_min=(lindemann_criteria, 0.0), out_file_name=f'{jiggle_path}/supercell_min.xyz')
                            shutil.copy(f'{parent}/energy.inp', jiggle_path)
                            extract_comment_change_inp(input_file = f'{jiggle_path}/energy.inp', xyz_file = f'{jiggle_path}/supercell_min.xyz')
                            
            print('Folder set-up complete!')
            print('Rerun script to begin simulations')


        else:  #case where files arent yet created and minmisation is not complete
            print('ERROR: minimisation step has not been performed, is it currently queued? You may need to manually run minimisation.')
            sys.exit()

    else: # Files are already created
        os.chdir(parent)
        progress = os.path.join(parent, 'progress.txt')
        if not os.path.exists(progress):
            print('Creating progress file...')
            #lists for volumes, strains, and jiggles

            strain_list = [n for n in range (1, (n_strain + 1))]

            jiggle_list = [n for n in range (1, (n_jiggle + 1))]

            column_widths = [15, 15, 15, 15, 15]

            total_width = sum(column_widths) + len(column_widths) * 4 + 1  #Total seperating line widths

            format_string = "| {:^{w0}} | {:^{w1}} | {:^{w2}} | {:^{w3}} |\n".format

            with open("progress.txt", "w") as f:
                f.write(format_string("", "", "","Complete?", w0=column_widths[0], w1=column_widths[1], w2=column_widths[2], w3=column_widths[3]))

                #Dividing line
                f.write(f"{'+' + '-' * (column_widths[0] + 2)}+{'-' * (column_widths[1] + 2)}+{'-' * (column_widths[2] + 2)}+{'-' * (column_widths[3] + 2)}+\n")

                for vol in vol_list:
                    f.write(format_string(f"vol{vol}", "", "", "", w0=column_widths[0], w1=column_widths[1], w2=column_widths[2], w3=column_widths[3]))
                    for strain in strain_list:
                        f.write(format_string("", f"strain{strain}", "", "", w0=column_widths[0], w1=column_widths[1], w2=column_widths[2], w3=column_widths[3]))
                        for jiggle in jiggle_list:
                            f.write(format_string("", "", f"jiggle{jiggle}", "n", w0=column_widths[0], w1=column_widths[1], w2=column_widths[2], w3=column_widths[3]))

                    #Dividing line for each vol block
                    f.write(f"{'+' + '-' * (column_widths[0] + 2)}+{'-' * (column_widths[1] + 2)}+{'-' * (column_widths[2] + 2)}+{'-' * (column_widths[3] + 2)}+\n")
            print('Progress file created')


    ############################################################# RUN ENERGY_FORCE CALCULATION FOR EACH DIFFERENT DEFORMATION #################################################################################

        def count_jobs():

            result = subprocess.run(['squeue', '--me', '--format="%j"'],
                                            capture_output = True, text = True)
            if result.returncode == 0:
                job_names = [line for line in result.stdout.splitlines() if line]
                job_count = len(job_names) - 1
            return job_count


        index = 0
        #While there are still simulations to be run (there are columns with n in)
        job_count = count_jobs()
        max_jobs = 5 #16 nodes used in sims and max allocation per user of 1024. 64 maximum jobs, 60 submitted leaves room for other submission.

        if job_count >= max_jobs:
            print('Max jobs reached, exiting.')
            sys.exit()

        while job_count <= max_jobs:
            job_count = count_jobs()
            current_vol, current_strain, current_jiggle, simulation_complete_check, index = run_simulation(index)
            if simulation_complete_check == False:
                if (current_vol is None and current_strain is None and current_jiggle is None):
                    print('All required jobs complete or running!')
                    break
                time.sleep(1)  #avoid rapid job submission, sleep for 1 second
        print('Job submission completed!')
        sys.exit()

