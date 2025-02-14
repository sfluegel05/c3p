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
        if atom.GetAtomicNum() == 6 and atom.GetTotalDegree() == 3: # Find carbons with 3 bonds
            neighbors = [neighbor.GetAtomicNum() for neighbor in atom.GetNeighbors()]
            if neighbors.count(6) == 1 and neighbors.count(8) == 1 and neighbors.count(1) == 2: # One C, one O, and two Hs
               o_neighbor = [neighbor for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() == 8][0]
               if o_neighbor.GetTotalDegree() == 1:
                terminal_ch2oh_count += 1
    
    if terminal_ch2oh_count != 2:
       return False, "Molecule does not have exactly two terminal CH2OH groups"


    # Check for CHOH groups in the chain
    for atom in atoms:
        if atom.GetAtomicNum() == 6 and atom.GetTotalDegree() == 3: # carbons connected to 3 other atoms
            neighbors = [neighbor.GetAtomicNum() for neighbor in atom.GetNeighbors()]
            if neighbors.count(6) == 2 and neighbors.count(8) == 1 and neighbors.count(1) == 1: # two carbon neighbours and 1 oxygen and hydrogen
                o_neighbor = [neighbor for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() == 8][0]
                if o_neighbor.GetTotalDegree() != 1: # check if the oxygen is part of an OH group
                    return False, "Molecule does not have CHOH in the chain"
            elif neighbors.count(6) != 1 and neighbors.count(8) != 1 and neighbors.count(1) != 2:
              return False, "Molecule does not have CHOH in the chain (found an invalid carbon neighbor)"
            
    # Verify that there are only carbons, hydrogens, and oxygens
    for atom in atoms:
       if atom.GetAtomicNum() not in [1, 6, 8]:
            return False, "Molecule contains atoms other than C, H, and O"

    return True, "Molecule is an acyclic polyol with the general formula HOCH2[CH(OH)]nCH2OH"