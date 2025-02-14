"""
Classifies: CHEBI:18254 ribonucleoside
"""
from rdkit import Chem

def is_ribonucleoside(smiles: str):
    """
    Determines if a molecule is a ribonucleoside based on its SMILES string.
    A ribonucleoside is a nucleoside where the sugar component is D-ribose, which is linked to a nucleobase.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ribonucleoside, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pattern for D-ribose sugar
    ribose_pattern = Chem.MolFromSmarts("O[C@@H]1[C@H](O)[C@@H](O)[C@H](CO)O1")

    # Patterns for common nucleobases with potential N-glycosidic link with ribose
    # Purine nucleobase linked via N9
    purine_linked_pattern = Chem.MolFromSmarts("n1c[nH]c2ncnc12")
    # Pyrimidine nucleobase linked via N1
    pyrimidine_linked_pattern = Chem.MolFromSmarts("n1c(=O)[nH]c(=O)c1")

    # Check for D-ribose component in the molecule
    if not mol.HasSubstructMatch(ribose_pattern):
        return False, "D-ribose component not found in the molecule"
    
    has_connected_nucleobase = False
    
    # Find the D-ribose substructure matches in the molecule
    ribose_matches = mol.GetSubstructMatches(ribose_pattern)
    
    for match in ribose_matches:
        ribose_atoms = set(match)
        for atom_idx in ribose_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                # Check if the neighboring atoms match the nucleobase pattern
                if not neighbor_idx in ribose_atoms:
                    if (mol.HasSubstructMatch(purine_linked_pattern, atomIdx=neighbor_idx) or
                        mol.HasSubstructMatch(pyrimidine_linked_pattern, atomIdx=neighbor_idx)):
                        has_connected_nucleobase = True
                        break
            if has_connected_nucleobase:
                break

    if not has_connected_nucleobase:
        return False, "Valid nucleobase is not connected to D-ribose"

    return True, "Valid ribonucleoside with D-ribose linked to a nucleobase."