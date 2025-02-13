"""
Classifies: CHEBI:36702 2-acyl-1-alkyl-sn-glycero-3-phosphocholine
"""
from rdkit import Chem

def is_2_acyl_1_alkyl_sn_glycero_3_phosphocholine(smiles: str):
    """
    Determines if a molecule is a 2-acyl-1-alkyl-sn-glycero-3-phosphocholine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-acyl-1-alkyl-sn-glycero-3-phosphocholine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pattern for phosphocholine group
    phosphocholine_pattern = Chem.MolFromSmarts("[O-]P(=O)(OC)OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "No phosphocholine group found"

    # Pattern for sn-1 alkyl ether linkage
    alkyl_chain_pattern = Chem.MolFromSmarts("COCCCC")
    ether_matches = mol.GetSubstructMatches(alkyl_chain_pattern)
    if len(ether_matches) < 1:
        return False, "No alkyl chain found (sn-1 position)"

    # Pattern for sn-2 acyl group in ester linkage
    acyl_group_pattern = Chem.MolFromSmarts("OC(=O)C")
    acyl_matches = mol.GetSubstructMatches(acyl_group_pattern)
    if len(acyl_matches) < 1:
        return False, "No acyl group found (sn-2 position)"

    return True, "Molecule matches the 2-acyl-1-alkyl-sn-glycero-3-phosphocholine structure"