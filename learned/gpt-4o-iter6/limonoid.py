"""
Classifies: CHEBI:39434 limonoid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

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
    
    # Check high oxygen content
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 6:  # Arbitrary threshold based on provided examples
        return False, "Insufficient oxygenation for typical limonoid"
    
    # Check for core triterpenoid structure (approximate by ring count and atom count)
    ring_count = Chem.GetSSSR(mol)
    if ring_count < 4:
        return False, "Insufficient number of rings for a triterpenoid structure"
    
    # Estimate presence of furan group (five-membered oxygen-containing aromatic ring)
    furan_pattern = Chem.MolFromSmarts("c1ccoc1")
    if not mol.HasSubstructMatch(furan_pattern):
        return False, "Lacks furan ring typical in limonoids"

    # Estimate molecular weight to ensure the size is typical for triterpenoids
    mol_wt = Descriptors.MolWt(mol)
    if mol_wt < 400 or mol_wt > 700:
        return False, "Molecular weight is atypical for limonoids"

    return True, "Structure consistent with limonoid features"