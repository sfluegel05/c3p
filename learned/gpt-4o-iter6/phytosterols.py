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

    # Sterane core: refined detection of the tetracyclic structure, specific for ring systems
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C(C1)CCC3C2CCC4=C3C=CCC4")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No tetracyclic steroid backbone detected"

    # Hydroxy group at specific position typically in the ring systems
    hydroxy_pattern = Chem.MolFromSmarts("[C@@H]1CC[C@@H]2[C@]1(CC[C@@]1([C@@]2(CC[C@]2([C@@]1(C)CCC=C2)O))C)C")
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "Typical hydroxyl group not found at specific position"

    # Check for side chain variations specifically in higher phytosterols
    side_chain_variations = [
        Chem.MolFromSmarts("CC(C)CCC=C"),  # typical for variations in side chain like extra methyl or double bond
    ]
    has_variations = any(mol.HasSubstructMatch(patt) for patt in side_chain_variations)
    
    # Condition of unsaturation or side chain changes which is typically for phytosterols
    if not has_variations:
        return False, "No variation in side chain length or unsaturation typical of phytosterols"

    return True, "Contains tetracyclic steroid backbone with phytosterol-specific side chain variations"