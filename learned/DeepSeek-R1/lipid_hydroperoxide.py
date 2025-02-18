"""
Classifies: CHEBI:61051 lipid hydroperoxide
"""
"""
Classifies: CHEBI:64606 lipid hydroperoxide
"""
from rdkit import Chem

def is_lipid_hydroperoxide(smiles: str):
    """
    Determines if a molecule is a lipid hydroperoxide based on its SMILES string.
    A lipid hydroperoxide is a lipid carrying one or more hydroperoxy (-OOH) groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lipid hydroperoxide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Check for hydroperoxy groups (-OOH)
    hydroperoxy_pattern = Chem.MolFromSmarts("[OX2]-[OX2H]")
    hp_matches = mol.GetSubstructMatches(hydroperoxy_pattern)
    if not hp_matches:
        return False, "No hydroperoxy (-OOH) groups found"
    
    # Check if it's a lipid: long carbon chain (at least 10 carbons) 
    # and possibly a carboxylic acid/ester (common in fatty acids)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 10:
        return False, f"Only {carbon_count} carbons, likely not a lipid"
    
    # Check for lipid-like features: carboxylic acid or ester groups
    carboxylic_acid = Chem.MolFromSmarts("C(=O)O")
    ester = Chem.MolFromSmarts("C(=O)O[C]")
    if not mol.HasSubstructMatch(carboxylic_acid) and not mol.HasSubstructMatch(ester):
        return False, "No carboxylic acid or ester group, not a typical lipid"
    
    return True, f"Lipid with {len(hp_matches)} hydroperoxy group(s)"