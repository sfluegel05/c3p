"""
Classifies: CHEBI:18254 ribonucleoside
"""
from rdkit import Chem

def is_ribonucleoside(smiles: str):
    """
    Determines if a molecule is a ribonucleoside based on its SMILES string.
    A ribonucleoside is a nucleoside where the sugar component is D-ribose.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ribonucleoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES to RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Substructure pattern for D-ribose 
    ribose_pattern = Chem.MolFromSmarts("O[C@@H]1[C@H](O)[C@@H](O)[C@H](CO)O1")

    if not mol.HasSubstructMatch(ribose_pattern):
        return False, "D-ribose component not found"

    # More detailed purine pattern including common modifications
    purine_pattern = Chem.MolFromSmarts("c1ncnc2[nH]cnc12")
    
    # More detailed pyrimidine pattern including common modifications
    pyrimidine_pattern = Chem.MolFromSmarts("c1c[nH]cnc1(=O)")
    
    # Ensure each matching pattern is connected directly to ribose
    substructs = mol.GetSubstructMatches(ribose_pattern)
    has_connected_nucleobase = False

    for substruct in substructs:
        ribose_idxs = set(substruct)  # Indices of the ribose substructure
        for atom in mol.GetAtoms():
            if atom.GetIdx() in ribose_idxs:
                for neighbor in atom.GetNeighbors():
                    connected_nucleobase = any(mol.HasSubstructMatch(pattern, atom.GetIdx())
                                               for pattern in [purine_pattern, pyrimidine_pattern])
                    if connected_nucleobase:
                        has_connected_nucleobase = True
                        break
            if has_connected_nucleobase:
                break

    if not has_connected_nucleobase:
        return False, "Nucleobase pattern not discernible or not connected to D-ribose"

    return True, "Molecule contains a D-ribose moiety linked to a nucleobase, classified as a ribonucleoside"