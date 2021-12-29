"""
Programmer: Aamir Alaud Din, PhD
File: lid_functions.py
Date: 2021.12.28

License:
    This module contains functions to obtain data from various
    sections of LAMMPS input data (Lid) file.
    Copyright (C) 2021 Aamir Alaud Din, PhD

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY: without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see <https://www.gnu.org/licenses/>
"""


from colorama import Fore


def lid_header(fname):
    """ This function extracts header from LAMMPS input data file.

    Syntax:
        mol_data, typ_data, bnd_data = lid_header(fname)

    Input argument(s):
        fname: Name of LAMMPS input data file.

    Output argument(s):
        molecular_data: No. of atoms, bonds, angles, dihedrals, and impropers
        in LAMMPS input data file.

        types_data: No. of atom, bond, angle, dihedral, and improper types in
        LAMMPS input data file.

        boundary_data: Data of periodic boundary box which contains lower and
        upper limits of simulation box in x, y, and z directions.

    Example:
        mol_data, typ_data, bnd_data = lid_header("system.data")

    """

    with open(fname, "r") as ifile:
        lines = ifile.readlines()
        ifile.close()

    sections = ["Masses", "Atoms", "Bonds", "Angles", "Dihedrals", "Impropers"]
    molecular_prop = ["atoms", "bonds", "angles", "dihedrals", "impropers"]
    types_prop = ["atom", "bond", "angle", "dihedral", "improper"]
    boundary_prop = ["xlo", "xhi", "ylo", "yhi", "zlo", "zhi"]
    molecular_data = [0, 0, 0, 0, 0]
    types_data = [0, 0, 0, 0, 0]
    boundary_data = [0, 0, 0, 0, 0, 0]

    for line in lines:
        line = line.split()

        if (len(line) >= 2 and line[1] in molecular_prop):
            molecular_data[molecular_prop.index(line[1])] = line[0]

        elif (len(line) >= 3 and line[1] in types_prop):
            types_data[types_prop.index(line[1])] = line[0]

        elif (len(line) >= 4 and line[2] in boundary_prop):
            boundary_data[boundary_prop.index(line[2])] = line[0]
            boundary_data[boundary_prop.index(line[2]) + 1] = line[1]

        elif ((len(line) >= 1) and (line[0] in sections)):
            break

    return molecular_data, types_data, boundary_data


def parse_lid(fname):
    """This function parses LAMMPS input data file and returns the line numbers
    at which various sections of LAMMPS input data start.

    Syntax:
        lid_sections = parse_lid(fname)

    Input argument(s):
        fname: Name of LAMMPS input data file.

    Output argument(s):
        The output argument of this function is a list containing six sublists
        i.e., 1. masses section, 2. atoms section, 3. bonds section,
        4. angles section, 5. dihedrals section, and 6. impropers section.
        Each sublist has two elements in which first element is the line
        number of LAMMPS input data file at which the corresponding section
        starts and second element is always zero. If first element of a sublist
        is zero, it means the corresponding section does not exist in LAMMPS
        input data file. Second zero element in any sublist is utilized by
        other functions of this module to extract required data of a section
        from LAMMPS input data file.

    Example:
        lid_sections = parse_lid("system.data")

    """

    with open(fname, "r") as ifile:
        lines = ifile.readlines()
        ifile.close()

    mas_section = [0, 0]
    atm_section = [0, 0]
    bnd_section = [0, 0]
    ang_section = [0, 0]
    dih_section = [0, 0]
    imp_section = [0, 0]

    counter = 0

    for line in lines:
        line = line.split()

        if (len(line) >= 1 and line[0] == "Masses"):
            mas_section[0] = counter

        elif (len(line) >= 1 and line[0] == "Atoms"):
            atm_section[0] = counter

        elif (len(line) >= 1 and line[0] == "Bonds"):
            bnd_section[0] = counter

        elif (len(line) >= 1 and line[0] == "Angles"):
            ang_section[0] = counter

        elif (len(line) >= 1 and line[0] == "Dihedrals"):
            dih_section[0] = counter

        elif (len(line) >= 1 and line[0] == "Impropers"):
            imp_section[0] = counter

        counter += 1

    return [mas_section, atm_section, bnd_section, ang_section, dih_section,
            imp_section]


def masses(fname):
    """This function returns masses section from LAMMPS input data file.

    Syntax:
        masses_data = masses(fname)

    Input argument(s):
        fname: Name of LAMMPS input data file.

    Output argument(s):
        masses_data: Data of atomic masses of elements in the relevant force
        field in LAMMPS input data file.

    Example:
        masses_data = masses("system.data")

    """

    in_data = parse_lid(fname)

    stop = []

    if in_data[0][0] == 0:
        print(Fore.RED + "Warning: No masses data available.")
        print(Fore.RED + "LAMMPS simulation can't run.")
        masses_data = None
        start = -1
        stop = -1

    elif in_data[0][0] == max(in_data)[0]:
        start = in_data[0][0] + 1
        stop = None
    else:
        for element in enumerate(in_data):

            if element[1][0] - in_data[0][0] == 0:
                stop.append(5000000000000)

            elif element[1][0] - in_data[0][0] < 0:
                stop.append(5000000000000)

            elif element[1][0] - in_data[0][0] > 0:
                stop.append(element[1][0] - in_data[0][0])

        stop = stop.index(min(stop))
        stop = in_data[stop][0]
        start = in_data[0][0] + 1

    with open(fname, "r") as infile:
        lines = infile.readlines()[start:stop]
        infile.close()

    masses_data = []
    for line in lines:
        line = line.split()

        if line != []:
            masses_data.append(line)

    return masses_data


def atoms(fname):
    """This function returns atoms section from LAMMPS input data file.

    Syntax:
        atoms_data = atoms(fname)

    Input argument(s):
        fname: Name of LAMMPS input data file.

    Output argument(s):
        atoms_data: Atoms data (atomic IDs, atom types, molecular IDs, charge,
        and x, y, z coordinates of all atoms in the simulation box) from
        LAMMPS input data file.

    Example:
        atoms_data = atoms("system.data")

    """

    in_data = parse_lid(fname)

    stop = []

    if in_data[1][0] == 0:
        print(Fore.RED + "Warning: No atoms data available.")
        print(Fore.RED + "LAMMPS simulation can't run.")
        atoms_data = None
        start = -1
        stop = -1

    elif in_data[1][0] == max(in_data)[0]:
        start = in_data[1][0] + 1
        stop = None

    else:
        for element in enumerate(in_data):

            if element[1][0] - in_data[1][0] == 0:
                stop.append(5000000000000)

            elif element[1][0] - in_data[1][0] < 0:
                stop.append(5000000000000)

            elif element[1][0] - in_data[1][0] > 0:
                stop.append(element[1][0] - in_data[1][0])

        stop = stop.index(min(stop))
        stop = in_data[stop][0]
        start = in_data[1][0] + 1

    with open(fname, "r") as infile:
        lines = infile.readlines()[start:stop]
        infile.close()

    atoms_data = []
    for line in lines:
        line = line.split()

        if line != []:
            atoms_data.append(line)

    return atoms_data


def bonds(fname):
    """This function returns bonds section from LAMMPS input data file.

    Syntax:
        bonds_data = bonds(fname)

    Input argument(s):
        fname: Name of LAMMPS input data file.

    Output argument(s):
        bonds_data: Bonds data (bond IDs, bond types, atomic IDs of bonding
        atoms) from LAMMPS input data file.

    Example:
        bonds_data = bonds("system.data")

    """

    in_data = parse_lid(fname)

    stop = []

    if in_data[2][0] == 0:
        print(Fore.RED + "Warning: No atoms data available.")
        print(Fore.RED + "LAMMPS simulation may produce wrong results.")
        bonds_data = None
        start = -1
        stop = -1

    elif in_data[2][0] == max(in_data)[0]:
        start = in_data[2][0] + 1
        stop = None

    else:
        for element in enumerate(in_data):

            if element[1][0] - in_data[2][0] == 0:
                stop.append(5000000000000)

            elif element[1][0] - in_data[2][0] < 0:
                stop.append(5000000000000)

            elif element[1][0] - in_data[2][0] > 0:
                stop.append(element[1][0] - in_data[2][0])

        stop = stop.index(min(stop))
        stop = in_data[stop][0]
        start = in_data[2][0] + 1

    with open(fname, "r") as infile:
        lines = infile.readlines()[start:stop]
        infile.close()

    bonds_data = []
    for line in lines:
        line = line.split()

        if line != []:
            bonds_data.append(line)

    return bonds_data


def angles(fname):
    """This function returns angles section from LAMMPS input data file.

    Syntax:
        angles_data = angles(fname)

    Input argument(s):
        fname: Name of LAMMPS input data file.

    Output argument(s):
        angles_data: Angles data (angle IDs, angle types, and atomic IDs having
        angle) from LAMMPS input data file.

    Example:
        angles_data = angles("system.data")

    """

    in_data = parse_lid(fname)

    stop = []

    if in_data[3][0] == 0:
        print(Fore.RED + "Warning: No angles data available.")
        print(Fore.RED + "LAMMPS simulation may produce wrong results.")
        angles_data = None
        start = -1
        stop = -1

    elif in_data[3][0] == max(in_data)[0]:
        start = in_data[3][0] + 1
        stop = None

    else:
        for element in enumerate(in_data):

            if element[1][0] - in_data[3][0] == 0:
                stop.append(5000000000000)

            elif element[1][0] - in_data[3][0] < 0:
                stop.append(5000000000000)

            elif element[1][0] - in_data[3][0] > 0:
                stop.append(element[1][0] - in_data[3][0])

        stop = stop.index(min(stop))
        stop = in_data[stop][0]
        start = in_data[3][0] + 1

    with open(fname, "r") as infile:
        lines = infile.readlines()[start:stop]
        infile.close()

    angles_data = []
    for line in lines:
        line = line.split()

        if line != []:
            angles_data.append(line)

    return angles_data


def dihedrals(fname):
    """This function returns dihedrals section from LAMMPS input data file.

    Syntax:
        dihedrals_data = dihedrals(fname)

    Input argument(s):
        fname: Name of LAMMPS input data file.

    Output argument(s):
        dihedrals_data: Dihedrals data (dihedral IDs, dihedral types, and
        atomic IDs having dihedral angle) from LAMMPS input data file.

    Example:
        dihedrals_data = dihedrals("system.data")

    """

    in_data = parse_lid(fname)

    stop = []

    if in_data[4][0] == 0:
        print(Fore.RED + "Warning: No dihedrals data available.")
        print(Fore.RED + "LAMMPS simulation may produce wrong results.")
        dihedrals_data = None
        start = -1
        stop = -1

    elif in_data[4][0] == max(in_data)[0]:
        start = in_data[4][0] + 1
        stop = None

    else:
        for element in enumerate(in_data):

            if element[1][0] - in_data[4][0] == 0:
                stop.append(5000000000000)

            elif element[1][0] - in_data[4][0] < 0:
                stop.append(5000000000000)

            elif element[1][0] - in_data[4][0] > 0:
                stop.append(element[1][0] - in_data[4][0])

        stop = stop.index(min(stop))
        stop = in_data[stop][0]
        start = in_data[4][0] + 1

    with open(fname, "r") as infile:
        lines = infile.readlines()[start:stop]
        infile.close()

    dihedrals_data = []
    for line in lines:
        line = line.split()

        if line != []:
            dihedrals_data.append(line)

    return dihedrals_data


def impropers(fname):
    """This function returns impropers section from LAMMPS input data file.

    Syntax:
        impropers_data = impropers(fname)

    Input argument(s):
        fname: Name of LAMMPS input data file.

    Output argument(s):
        impropers_data: Impropers data (improper IDs, improper types, and
        atomic IDs having improper angle) from LAMMPS input data file.

    Example:
        impropers_data = impropers("system.data")

    """

    in_data = parse_lid(fname)

    stop = []

    if in_data[5][0] == 0:
        print(Fore.RED + "Warning: No impropers data available.")
        print(Fore.RED + "LAMMPS simulation may produce wrong results.")
        impropers_data = None
        start = -1
        stop = -1

    elif in_data[5][0] == max(in_data)[0]:
        start = in_data[5][0] + 1
        stop = None

    else:
        for element in enumerate(in_data):

            if element[1][0] - in_data[5][0] == 0:
                stop.append(5000000000000)

            elif element[1][0] - in_data[5][0] < 0:
                stop.append(5000000000000)

            elif element[1][0] - in_data[5][0] > 0:
                stop.append(element[1][0] - in_data[5][0])

        stop = stop.index(min(stop))
        stop = in_data[stop][0]
        start = in_data[5][0] + 1

    with open(fname, "r") as infile:
        lines = infile.readlines()[start:stop]
        infile.close()

    impropers_data = []
    for line in lines:
        line = line.split()

        if line != []:
            impropers_data.append(line)

    return impropers_data
