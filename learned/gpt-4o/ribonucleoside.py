"""
Classifies: CHEBI:18254 ribonucleoside
"""
from rdkit import Chem

def is_ribonucleoside(smiles: str):
    """
    Determines if a molecule is a ribonucleoside based on its SMILES string.
    A ribonucleoside is a nucleoside where the sugar component is D-ribose, 
    linked to a nucleobase via an N-glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ribonucleoside, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pattern for D-ribose sugar with a linkage capable position
    ribose_pattern = Chem.MolFromSmarts("OC[C@H]1O[C@H](CO)[C@@H](O)[C@H]1O")

    # Patterns for nucleobases
    purine_nucleobase = Chem.MolFromSmarts("c1[nH]c2ncnc2c1")  # Simplified purine core with linkage
    pyrimidine_nucleobase = Chem.MolFromSmarts("c1[nH]cnc(=O)c1")  # Simplified pyrimidine core with linkage

    # Check for D-ribose component in the molecule
    if not mol.HasSubstructMatch(ribose_pattern):
        return False, "D-ribose component not found in the molecule"
    
    # Check for nucleobase component and ensure linkage
    ribose_matches = mol.GetSubstructMatches(ribose_pattern)
    for ribose_match in ribose_matches:
        ribose_atoms = set(ribose_match)
        for atom_idx in ribose_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                # Ensure the neighbor is not part of ribose
                if neighbor_idx not in ribose_atoms:
                    # Check purine linkage
                    if mol.HasSubstructMatch(purine_nucleobase, atomIdx=neighbor_idx):
                        return True, "Valid ribonucleoside with purine base linked to D-ribose"
                    # Check pyrimidine linkage
                    if mol.HasSubstructMatch(pyrimidine_nucleobase, atomIdx=neighbor_idx):
                        return True, "Valid ribonucleoside with pyrimidine base linked to D-ribose"

    return False, "Valid nucleobase is not connected to D-ribose"