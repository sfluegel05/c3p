"""
Classifies: CHEBI:26125 phytosterols
"""
from rdkit import Chem

def is_phytosterols(smiles: str):
    """
    Determines if a molecule is a phytosterol based on its SMILES string.
    Phytosterols are sterols similar to cholesterol which occur in plants 
    and may vary in carbon side chains and/or presence or absence of a double bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phytosterol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generalized tetracyclic steroid core pattern for phytosterols
    steroid_core_pattern = Chem.MolFromSmarts("C1CC[C@H]2[C@@H]1CC[C@@H]3[C@H]2CC[C@H]4=C3[CH2]CC=C4")
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "No tetracyclic steroid backbone detected"

    # Hydroxyl groups typically present on sterol structures
    hydroxy_group_pattern = Chem.MolFromSmarts("O[C@H]([C@H]1CC[C@H]2[C@@H]1CCC3[C@H]2CCC4=C3[CH2]CC=C4)C")
    if not mol.HasSubstructMatch(hydroxy_group_pattern):
        return False, "Typical hydroxyl group on sterol backbone not detected"

    # Side chain variability: broader pattern allowing for different side chain lengths and unsaturations
    side_chain_variation_patterns = [
        Chem.MolFromSmarts("CC(C)C"),  # side chains with extra branches like methyl groups
        Chem.MolFromSmarts("C=C"),     # side chains with double bonds
    ]
    
    # Check for presence of side chain modifications
    has_side_chain_variation = any(mol.HasSubstructMatch(pattern) for pattern in side_chain_variation_patterns)
    if not has_side_chain_variation:
        return False, "No variation in side chain typical of phytosterols"

    return True, "Contains tetracyclic steroid backbone with phytosterol-specific side chain variations"