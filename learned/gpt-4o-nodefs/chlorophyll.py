"""
Classifies: CHEBI:28966 chlorophyll
"""
from rdkit import Chem

def is_chlorophyll(smiles: str):
    """
    Determines if a molecule is chlorophyll based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is chlorophyll, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for magnesium (Mg) ion coordinated in the molecule
    mg_atom = any(atom.GetAtomicNum() == 12 for atom in mol.GetAtoms())
    if not mg_atom:
        return False, "No magnesium atom found"

    # Recognize characteristics of a porphyrin ring with magnesium
    porphyrin_pattern = Chem.MolFromSmarts("[n]1[c][n]2[c][n]3[c][n]4[c][N+]1([Mg--])[n]5[c][n]6[c][n]7[c][N+]3([Mg--])[n]42")
    if not mol.HasSubstructMatch(porphyrin_pattern):
        return False, "No porphyrin-like structure found"

    # Check for long hydrocarbon chains typically associated with chlorophyll molecules
    hydrocarbon_chain_pattern = Chem.MolFromSmarts("CCCCCCCCCCCC")
    if not mol.HasSubstructMatch(hydrocarbon_chain_pattern):
        return False, "Long hydrocarbon chains characteristic of chlorophyll not found"
    
    # Check for ester groups that are typically found in chlorophyll
    ester_pattern = Chem.MolFromSmarts("COC(=O)")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 1:
        return False, "No ester groups found, common in chlorophyll"

    return True, "Matches chlorophyll characteristics: porphyrin-like ring with magnesium, hydrocarbon chains, and ester groups"