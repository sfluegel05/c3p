"""
Classifies: CHEBI:18000 aralkylamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aralkylamine(smiles: str):
    """
    Determines if a molecule is an aralkylamine based on its SMILES string.
    An aralkylamine is an alkylamine where at least one of the alkyl substituents is an aromatic group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aralkylamine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check for a basic amine (primary, secondary, or tertiary, but not quaternary)
    amine_pattern = Chem.MolFromSmarts("[N;H0,H1,H2]") # N with 0, 1 or 2 H
    if not mol.HasSubstructMatch(amine_pattern):
        return False, "No amine group found"

    # 2. Check for an aromatic ring (specifically, benzene ring)
    aromatic_ring_pattern = Chem.MolFromSmarts("c1ccccc1")
    if not mol.HasSubstructMatch(aromatic_ring_pattern):
        return False, "No aromatic ring found"

    # 3. Check for an alkyl bridge (chain of at least 1 carbon between amine and aromatic)
    # We need to make sure there's a chain of carbons between the N and aromatic carbons
    alkyl_bridge_pattern1 = Chem.MolFromSmarts("[N;H0,H1,H2]-[CX4]-[c]")
    alkyl_bridge_pattern2 = Chem.MolFromSmarts("[N;H0,H1,H2]-[CX4]-[CX4]-[c]")
    alkyl_bridge_pattern3 = Chem.MolFromSmarts("[N;H0,H1,H2]-[CX4]-[CX4]-[CX4]-[c]")
    # We can check for bridges from length 1 to 3, might be necessary to check for more

    if not (mol.HasSubstructMatch(alkyl_bridge_pattern1) or mol.HasSubstructMatch(alkyl_bridge_pattern2) or mol.HasSubstructMatch(alkyl_bridge_pattern3)):
        return False, "No alkyl chain connecting amine and aromatic ring"
        
    
    # If all conditions are met, then it is an aralkylamine
    return True, "Contains an amine group connected to an aromatic ring via an alkyl chain"