"""
Classifies: CHEBI:26125 phytosterols
"""
"""
Classifies: CHEBI:???? phytosterols
"""
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

def is_phytosterols(smiles: str):
    """
    Determines if a molecule is a phytosterol based on its SMILES string.
    Phytosterols are sterols similar to cholesterol with variations in the side chain and/or double bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phytosterol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Check for hydroxyl group on a ring atom (C3 position in steroid nucleus)
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]-[C@]")  # Approximate, may need adjustment
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl group on a ring atom"
    
    # Check for four fused rings (steroid nucleus)
    # Using a SMARTS pattern for the tetracyclic system (simplified)
    steroid_core = Chem.MolFromSmarts("[C]1[C@@]2([C@H]([C@H]3[C@@]([C@@]([C@]3(CC2)C)(C)CC1)C)CC4)C[C@@H]4O")
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid nucleus detected"
    
    # Check for side chain with at least 8 carbons and either a double bond or ethyl group
    side_chain_double_bond = Chem.MolFromSmarts("C=C")
    ethyl_group = Chem.MolFromSmarts("[CH2CH3]")
    
    has_double = mol.HasSubstructMatch(side_chain_double_bond)
    has_ethyl = mol.HasSubstructMatch(ethyl_group)
    
    if not (has_double or has_ethyl):
        return False, "Side chain lacks double bond or ethyl group"
    
    # Molecular weight check (phytosterols are typically >370 Da)
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt < 370:
        return False, f"Molecular weight too low ({mol_wt:.1f} < 370)"
    
    return True, "Steroid nucleus with hydroxyl group and plant-specific side chain"