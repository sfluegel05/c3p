"""
Classifies: CHEBI:39434 limonoid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdqueries

def is_limonoid(smiles: str):
    """
    Determines if a molecule is a limonoid based on its SMILES string.
    A limonoid is a highly oxygenated triterpenoid with a 4,4,8-trimethyl-17-furanylsteroid structure.

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
    
    # Check high oxygen content (more refined threshold based on example data)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 10:  # Increase threshold; must review example spectrum
        return False, "Insufficient oxygenation for typical limonoid"
    
    # Attempt to identify core skeletal features (e.g., highly fused ring system)
    # Placeholders for now, as SMARTS-specific motifs would need derivation
    steroid_skeleton = Chem.MolFromSmarts("C1CCC2C1()C3C2CCC4C3=CC5C4C(O)=CC6C5(CC=C6)=O")
    if not mol.HasSubstructMatch(steroid_skeleton):
        return False, "Lacks steroid-type core"

    # Identify presence of furan ring
    furan_pattern = Chem.MolFromSmarts("c1ccoc1")
    if not mol.HasSubstructMatch(furan_pattern):
        return False, "Lacks furan ring typical in limonoids"
    
    # Molecular weight for triterpenoids, reflecting highly oxygenated nature
    mol_wt = Descriptors.MolWt(mol)
    if mol_wt < 450 or mol_wt > 1200:  # Broad range based on structural diversity
        return False, "Molecular weight is atypical for limonoids"

    return True, "Structure consistent with limonoid features"