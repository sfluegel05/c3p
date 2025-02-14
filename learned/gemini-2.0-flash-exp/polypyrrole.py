"""
Classifies: CHEBI:38077 polypyrrole
"""
from rdkit import Chem
from rdkit.Chem import Draw

def is_polypyrrole(smiles: str):
    """
    Determines if a molecule is a polypyrrole based on its SMILES string.
    A polypyrrole is a compound composed of two or more pyrrole units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polypyrrole, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the pyrrole substructure using SMARTS (more specific)
    pyrrole_smarts = "[nX1H1]1cc[cX3]c1"
    pyrrole_pattern = Chem.MolFromSmarts(pyrrole_smarts)

    # Get all matches of the pattern
    pyrrole_matches = mol.GetSubstructMatches(pyrrole_pattern)
    
    if len(pyrrole_matches) == 0:
        return False, "Contains no pyrrole units, therefore not a polypyrrole."

    # Count the number of pyrrole rings.
    num_pyrroles = 0

    #Initialize a set to track all nitrogen atoms in all pyrrole rings
    pyrrole_n_atoms = set()

    for match in pyrrole_matches:
        #Check if a nitrogen atom from the current match is in the already matched set of N atoms.
        nitrogen_atom_index = -1
        for atom_index in match:
            atom = mol.GetAtomWithIdx(atom_index)
            if atom.GetAtomicNum() == 7:
                nitrogen_atom_index = atom_index
                break #Only one N atom per pyrrole ring
        
        if nitrogen_atom_index != -1 and nitrogen_atom_index not in pyrrole_n_atoms:
            num_pyrroles += 1
            pyrrole_n_atoms.add(nitrogen_atom_index)


    # Classify based on the count of pyrrole rings
    if num_pyrroles >= 2:
        return True, f"Contains {num_pyrroles} pyrrole units, therefore a polypyrrole."
    else:
        return False, f"Contains only {num_pyrroles} pyrrole unit(s), therefore not a polypyrrole."