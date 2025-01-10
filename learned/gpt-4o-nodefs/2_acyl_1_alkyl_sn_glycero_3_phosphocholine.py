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

    # Revised pattern for the phosphocholine group with varying bonds
    phosphocholine_pattern = Chem.MolFromSmarts("[O-]P(=O)(OCC[N+](C)(C)C)O")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "No phosphocholine group found"

    # Updated flexible pattern for the ether linkage at sn-1 position with optional variations
    ether_chain_pattern = Chem.MolFromSmarts("CO[C@H](CO)([C,CX4])")
    if not mol.HasSubstructMatch(ether_chain_pattern):
        return False, "No ether linkage with alkyl chain found (sn-1 position)"

    # Flexible pattern for sn-2 ester-acyl group
    acyl_group_pattern = Chem.MolFromSmarts("OC(=O)[C,CX4]")
    if not mol.HasSubstructMatch(acyl_group_pattern):
        return False, "No acyl ester group found (sn-2 position)"

    return True, "Molecule matches the 2-acyl-1-alkyl-sn-glycero-3-phosphocholine structure"