"""
Classifies: CHEBI:25409 monoterpenoid
"""
"""
Classifies: CHEBI:36581 monoterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoterpenoid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid based on its SMILES string.
    Monoterpenoids are derived from monoterpenes (C10 skeleton) and may have structural modifications.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoterpenoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Check terpene classification using RDKit's terpene fingerprint
    try:
        fp = rdMolDescriptors.GetTerpeneFingerprint(mol)
    except:
        return False, "Error in terpene fingerprint calculation"
    
    # Check for monoterpene (m) or monoterpenoid (if present in fingerprint)
    if 'm' not in fp:
        return False, "Not classified as a monoterpene-derived structure"
    
    # Check carbon count considering possible modifications (original C10 ± some)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (8 <= c_count <= 20):
        return False, f"Carbon count {c_count} outside typical monoterpenoid range (8-20)"
    
    # Additional check for terpenoid nature (presence of oxygen is common but not mandatory)
    # Skipping as some monoterpenoids may be hydrocarbons
    
    return True, "Classified as monoterpenoid based on terpene fingerprint and carbon count"