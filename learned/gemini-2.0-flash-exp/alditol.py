"""
Classifies: CHEBI:17522 alditol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_alditol(smiles: str):
    """
    Determines if a molecule is an alditol based on its SMILES string.
    An alditol is an acyclic polyol with the general formula HOCH2[CH(OH)]nCH2OH.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alditol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for cycles
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings, not an alditol"

    # Get all atoms in the molecule
    atoms = mol.GetAtoms()

    # Check for the two terminal CH2OH groups
    terminal_ch2oh_count = 0
    for atom in atoms:
        if atom.GetAtomicNum() == 6 and atom.GetTotalDegree() == 3 : # carbons with three bonds
            neighbors = [neighbor for neighbor in atom.GetNeighbors()] # get neighbours
            oxygen_neighbor = [neigh for neigh in neighbors if neigh.GetAtomicNum() == 8] # get all O neighbours
            if len(oxygen_neighbor) == 1: # if there is one oxygen
                o_atom = oxygen_neighbor[0]
                if o_atom.GetTotalDegree() == 1: #oxygen needs to have a single bond
                    h_neighbors = [h for h in o_atom.GetNeighbors() if h.GetAtomicNum() == 1]
                    if len(h_neighbors) == 1: #oxygen must be bound to one hydrogen
                            h_neighbors = [neigh for neigh in neighbors if neigh.GetAtomicNum() == 1] #Get H neighbours from carbon
                            if len(h_neighbors) == 2: #There must be 2 H's on carbon
                                terminal_ch2oh_count += 1

    if terminal_ch2oh_count != 2:
        return False, f"Molecule does not have exactly two terminal CH2OH groups, found {terminal_ch2oh_count}"


    # Check for CHOH groups in the chain
    for atom in atoms:
        if atom.GetAtomicNum() == 6 and atom.GetTotalDegree() == 3: # carbons with three bonds
            neighbors = [neighbor for neighbor in atom.GetNeighbors()]
            oxygen_neighbor = [neigh for neigh in neighbors if neigh.GetAtomicNum() == 8] # get all O neighbours
            if len(oxygen_neighbor) == 1: # if there is one oxygen
                o_atom = oxygen_neighbor[0]
                if o_atom.GetTotalDegree() == 1: #oxygen needs to have a single bond
                    h_neighbors = [h for h in o_atom.GetNeighbors() if h.GetAtomicNum() == 1]
                    if len(h_neighbors) == 1: #oxygen must be bound to one hydrogen
                        
                        c_count = 0
                        for neigh in neighbors:
                            if neigh.GetAtomicNum() == 6:
                                c_count +=1
                        if c_count != 2:
                           return False, f"Molecule does not have CHOH in the chain, found carbon chain with carbon neighbours {c_count}"

    # Verify that there are only carbons, hydrogens, and oxygens
    for atom in atoms:
       if atom.GetAtomicNum() not in [1, 6, 8]:
            return False, "Molecule contains atoms other than C, H, and O"

    return True, "Molecule is an acyclic polyol with the general formula HOCH2[CH(OH)]nCH2OH"