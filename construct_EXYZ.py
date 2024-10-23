import numpy as np 
import os
import sys

class Construct_EXYZ:

    def __init__(self, xyz_file=None, force_file=None, output_file = None):

        files = [file for file in os.listdir()]

        # If xyz_file, output_file or force_file are not passed, find them in the directory
        if output_file is None:
            for file in files:
                if '.out' in file:
                    output_file = file
                    break

        if force_file is None:
            for file in files:
                if '_0.xyz' in file:
                    force_file = file
                    break
        
        if xyz_file is None:
            for file in files:
                if 'min.xyz' in file:
                    xyz_file = file
                    break

        # Error handling if files not found
        if force_file is None:
            print('ERROR: force file not found, ensure it contains the string "frc", or pass name manually (force_file =)')
            sys.exit()
        if xyz_file is None:
            print('ERROR: position file not found, ensure it contains the string "pos", or pass file name manually (xyz_file =)')
            sys.exit()
        if output_file is None:
            print('ERROR: cp2k output file not found, ensure it contains the string ".out", or pass name manually (output_file =)')
            sys.exit()

        # Set class attributes
        self.xyz_file = xyz_file
        self.force_file = force_file
        self.output_file = 'result.out'

        # Extract number of atoms from xyz_file
        with open(self.xyz_file, 'r') as file:
            lines = file.readlines()

            self.n_atoms = int(lines[0].split()[0])

    



    def extract_xyz(self, input_file = None):

        """
        Extracts data for the given steps. If steps are not provided, all steps are extracted.
        Parameters:
        - steps: Can be a list of specific steps, a single step (int), or a range of steps (tuple)
        - If file is not from MD (i.e. does not contain multiple steps), just extracts the xyz.
        """
        if input_file == None:
            input_file = self.xyz_file

        atoms = []
        positions = []
        with open(input_file, 'r') as file:
            lines = file.readlines()
        for line in lines[2:]:
            elements = line.split()
            atom, x, y, z = str(elements[0]), float(elements[1]), float(elements[2]), float(elements[3])
            atoms.append(atom)
            positions.append((x, y, z))

        lattice = lines[1].split(';')[2].split(':')[1].strip().replace('[', '').replace(']', '').replace(',', ' ')

        
        return atoms, positions, lattice


    
    def extract_forces(self, input_file = None):
        forces = []
        if input_file == None:
            input_file = self.force_file
        with open(input_file, 'r') as file:
            lines = file.readlines()
        non_empty_lines = [line for line in lines if line.strip()]
        for line in non_empty_lines[2:(2+self.n_atoms)]:
            elements = line.split()
            x_frc, y_frc, z_frc = float(elements[3]) * 51.422, float(elements[4]) * 51.422, float(elements[5]) * 51.422
            forces.append((x_frc, y_frc, z_frc))
        
      
        return forces
    
    def extract_energy(self, input_file = None):

        if input_file == None:
            input_file = self.output_file
        energy = None
        with open (input_file, 'r') as f:
            lines = f.readlines()
        for line in lines:
            if ' ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.]:' in line:
                energy = line.split()[8]
                break
        return energy


    def get_pos_force_ener(self):

        atoms, positions, lattice = self.extract_xyz()
        forces = self.extract_forces()
        energy = self.extract_energy()

        return positions, forces, energy, atoms, lattice
    
    
    def create_database(self):
        
        data = []
        

        positions, forces, energy, atoms, lattice = self.get_pos_force_ener()

        data.append(self.n_atoms)
        comment_string = f'dft_energy={energy} Lattice="{lattice}" Properties=species:S:1:pos:R:3:dft_force:R:3'
        data.append(comment_string)
        for i in range(self.n_atoms):
            string = f'{atoms[i]} {" ".join(map(str, positions[i]))} {" ".join(map(str, forces[i]))}'
            data.append(string)
        
        return data 
    


      
                






            







# 0.884675671002E+05 volume bohr**3, 0.4388419E+02   0.4493703E+02   0.4486128E+02 cell lengths bohr 
        

            

        

