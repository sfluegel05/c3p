"""
Classifies: CHEBI:18035 diglyceride
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_diglyceride(smiles: str):
    """
    Determines if a molecule is a diglyceride based on its SMILES string.
    A diglyceride is a glycerol molecule with two of its hydroxy groups acylated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a diglyceride, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Identify the three-carbon glycerol backbone pattern
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Identify the ester linkage pattern (C(=O)O)
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    # Check that there are exactly two ester linkages
    if len(ester_matches) != 2:
        return False, f"Found {len(ester_matches)} ester linkages, require exactly 2 for diglyceride"

    # Check attachment of ester linkages to the glycerol backbone
    # The expectation is that these ester linkages are attached to the primary carbon atoms of the glycerol backbone.
    match_found = False
    for ester_match in ester_matches:
        for atom_idx in ester_match:
            atom = mol.GetAtomWithIdx(atom_idx)
            # Check if the ester is attached to the glycerol backbone
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() in {atom.GetIdx() for m in mol.GetSubstructMatches(glycerol_pattern)}:
                    match_found = True
                    break
        if not match_found:
            return False, "Ester not bonded to glycerol backbone"

    # Check if there are two available positions for acyl groups, indicating a diglyceride
    if mol.GetSubstructMatches(glycerol_pattern).count(len(ester_matches)) != 2:
        return False, "More or less than two available acylation sites, not a diglyceride"

    return True, "Contains glycerol backbone with two fatty acid chains attached via ester bonds"