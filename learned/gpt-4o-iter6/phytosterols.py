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

    # Steroid backbone: detection of the sterane core structure
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C(C1)CCC3C2CCC4C3(CCCC4)C")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No tetracyclic steroid backbone detected"

    # Detect hydroxy group in specific positions
    hydroxy_pattern = Chem.MolFromSmarts("C[C@@H](O)")
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "Typical hydroxyl group not found"

    # Check for unsaturation: presence of double bonds in side chains or rings
    olefinic_patterns = [
        Chem.MolFromSmarts("C=C"),  # Broad olefinic bond pattern
    ]
    has_unsaturation = any(mol.HasSubstructMatch(patt) for patt in olefinic_patterns)
    
    # Typical side chain variations: noting longer chains or sp2 configurations 
    side_chain_variation = Chem.MolFromSmarts("C[C@H](C)CCC")
    if not mol.HasSubstructMatch(side_chain_variation) and not has_unsaturation:
        return False, "No variation in side chain length or unsaturation typical of phytosterols"

    return True, "Contains tetracyclic steroid backbone with phytosterol-specific side chain variations"