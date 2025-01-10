"""
Classifies: CHEBI:24043 flavones
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_flavones(smiles: str) -> (bool, str):
    """
    Determines if a molecule is a flavone based on its SMILES string.
    A flavone has a 2-aryl-1-benzopyran-4-one skeleton and its substituted derivatives.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a flavone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Flavone core structure - a 2-aryl-1-benzopyran-4-one
    flavone_pattern = Chem.MolFromSmarts("c1cc(ccc1)-c2c(oc3ccccc3c2=O)")[1]
    if not mol.HasSubstructMatch(flavone_pattern):
        return False, "No flavone core structure found"
    
    # Additional checks for substitutions can be performed here
    # Counting rings, aromatic systems can help check for possible derivatives
    try:
        num_aromatic_rings = Chem.Lipinski.NumAromaticRings(mol)
        if num_aromatic_rings < 2:
            return False, "Insufficient aromatic ring count for typical flavone derivatives"
    except:
        return None, None  # For difficult/corner cases we don't handle

    return True, "Contains flavone core structure with appropriate derivatives"