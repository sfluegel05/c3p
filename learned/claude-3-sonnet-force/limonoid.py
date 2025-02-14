"""
Classifies: CHEBI:39434 limonoid
"""
"""
Classifies: CHEBI:39572 limonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_limonoid(smiles: str):
    """
    Determines if a molecule is a limonoid based on its SMILES string.
    A limonoid is a highly oxygenated natural product derived from the tetracyclic triterpenoid skeleton,
    often containing a furan ring and specific methyl groups or other characteristic substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a limonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for furan ring
    furan_pattern = Chem.MolFromSmarts("c1ccoc1")
    has_furan = mol.HasSubstructMatch(furan_pattern)
    
    # Check for high oxygenation (at least 6 oxygens)
    n_oxygens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if n_oxygens < 6:
        return False, "Not enough oxygens for limonoid (need at least 6)"
    
    # Check for tetracyclic or rearranged tetracyclic skeleton
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    if n_rings < 4:
        return False, "Too few rings for limonoid (need at least 4)"
    
    # Check molecular weight - limonoids typically >400 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for limonoid"
    
    # Look for specific methyl groups or other characteristic substituents
    methyl_pattern = Chem.MolFromSmarts("[C]1(CCC2CCCC3CCCC(C2)C3C)C")
    has_methyl_groups = mol.HasSubstructMatch(methyl_pattern)
    
    # Classify based on the presence of key features
    if has_furan and has_methyl_groups:
        return True, "Contains furan ring, tetracyclic skeleton, high oxygenation, and characteristic methyl groups"
    elif has_furan:
        return True, "Contains furan ring, tetracyclic skeleton, and high oxygenation"
    elif has_methyl_groups:
        return True, "Contains tetracyclic skeleton, high oxygenation, and characteristic methyl groups"
    else:
        return True, "Contains tetracyclic skeleton and high oxygenation consistent with limonoid structure"