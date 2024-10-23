from pymatgen.ext.matproj import MPRester
from pymatgen.core.structure import Structure
from pymatgen.io.xyz import XYZ
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

api_key = ''  # Add api key here 
mpid = input('Input MP ID:')  #Input MP ID


with MPRester(api_key) as m:
    # Fetch the structure object
    struct = m.get_structure_by_material_id(f'mp-{mpid}')
    analyzer = SpacegroupAnalyzer(struct)
    struct = analyzer.get_conventional_standard_structure() #Retrun structure object as conventional cell (primitive is deafult)
    xyz_string = str(XYZ(struct))
    lattice = struct.lattice
    print(lattice)
    a, b, c = lattice.abc
    alpha, beta, gamma = lattice.angles
    lattice2 = struct.lattice.matrix
    vector_string = f'[[{lattice2[0][0]},{lattice2[0][1]},{lattice2[0][2]}],[{lattice2[1][0]},{lattice2[1][1]},{lattice2[1][2]}],[{lattice2[2][0]},{lattice2[2][1]},{lattice2[2][2]}]]'

    def classify_crystal_system(a, b, c, alpha, beta, gamma):
        if a != b != c and alpha != 90 and beta != 90 and gamma != 90:
            return "TRICLINIC"
        elif a != b != c and alpha == gamma == 90 and beta != 90:
            if gamma == 90 and alpha == 90:
                return "MONOCLINIC_GAMMA_AB"
            return "MONOCLINIC"
        elif a != b != c and alpha == beta == gamma == 90:
            return "ORTHORHOMBIC"
        elif a == b != c and alpha == beta == gamma == 90:
            if a == b and alpha == 90 and beta == 90 and gamma == 90:
                return "TETRAGONAL_AB"
            elif a == c and alpha == 90 and beta == 90 and gamma == 90:
                return "TETRAGONAL_AC"
            elif b == c and alpha == 90 and beta == 90 and gamma == 90:
                return "TETRAGONAL_BC"
            return "TETRAGONAL"
        elif a == b == c and alpha == beta == gamma != 90:
            return "RHOMBOHEDRAL"
        elif a == b != c and alpha == beta == 90 and (gamma == 60 or gamma == 120):
            if gamma == 60:
                return "HEXAGONAL_GAMMA_60"
            elif gamma == 120:
                return "HEXAGONAL_GAMMA_120"
            return "HEXAGONAL"
        elif a == b == c and alpha == beta == gamma == 90:
            return "CUBIC"
        else:
            return "UNKNOWN"

    crystal = classify_crystal_system(a, b, c, alpha, beta, gamma)

    comment = f"Lattice: a={a:.6f} b={b:.6f} c={c:.6f}; Angles: alpha={alpha:.6f} beta={beta:.6f} gamma={gamma:.6f}; Symmetry: {crystal}; Vectors: {vector_string}"
    xyz_lines = xyz_string.split('\n')
    xyz_lines.pop(1) #Remove original comment line
    xyz_lines.insert(1, comment)
    modified_xyz_string = '\n'.join(xyz_lines)


    if crystal == "UNKNOWN":
        print('WARNING, SYMMETRY IS NOT DETECTED')



    #write to XYZ file
    filename = "inp.xyz"
    with open(filename, "w") as f:
        f.write(modified_xyz_string)


