"""
Classifies: CHEBI:32955 epoxide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_epoxide(smiles: str):
    """
    Determines if a molecule is an epoxide based on its SMILES string.
    An epoxide is a three-membered ring containing an oxygen atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an epoxide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Define SMARTS pattern for epoxide (three-membered ring with one oxygen)
    epoxide_pattern = Chem.MolFromSmarts("C1OC1")

    # Check for a match of the pattern
    if mol.HasSubstructMatch(epoxide_pattern):
        # Check that the match is indeed a 3 member ring.
        matches = mol.GetSubstructMatches(epoxide_pattern)
        for match in matches:
            ring_atoms = [mol.GetAtomWithIdx(idx) for idx in match]
            ring_size = len(ring_atoms)
            if ring_size != 3:
                continue
            # check for ring bond count
            bonds_in_ring = 0
            for atom in ring_atoms:
                for neighbor in atom.GetNeighbors():
                    if neighbor in ring_atoms:
                        bonds_in_ring +=1
            if bonds_in_ring == 6: # total number of bonds in a 3 member ring
                return True, "Contains a three-membered ring with an oxygen"
            else:
                continue
        
        return False, "Does not contain an epoxide ring"

    return False, "Does not contain an epoxide ring"