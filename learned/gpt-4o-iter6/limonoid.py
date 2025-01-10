"""
Classifies: CHEBI:39434 limonoid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_limonoid(smiles: str):
    """
    Determines if a molecule is a limonoid based on its SMILES string.
    Limonoids are highly oxygenated triterpenoids with a specific core structure 
    featuring a 4,4,8-trimethyl-17-furanylsteroid skeleton and furan ring.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a limonoid, False otherwise
        str: Reason for classification
    """

    # Attempt to parse the SMILES string into a molecular structure
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of at least one furan ring (highly oxygenated domain)
    furan_pattern = Chem.MolFromSmarts("c1ccoc1")
    if not mol.HasSubstructMatch(furan_pattern):
        return False, "No furan ring identified, which is typical for limonoid"

    # Check for a pattern representing a highly oxygenated triterpenoid structure
    # Highly simplified suggestion that might not represent exact core, more definition needed
    core_pattern = Chem.MolFromSmarts("C1CC2(C)CC3C4(CCC5(C(OC)=O)C(O)=C)[C@]123C(=O)OC4=O")  # Example pattern
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Lacks core structural features typical of limonoids"

    # Advanced check for oxygenation extent, typically higher in oxygen atoms
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 8:  # Based on empirical examples that show extensive oxygenation
        return False, f"Insufficient oxygen atoms ({oxygen_count}) for high oxygenation typical of limonoids"

    mol_wt = Descriptors.MolWt(mol)
    if mol_wt < 400 or mol_wt > 1200:
        return False, "Molecular weight falls outside typical range for limonoids"

    return True, "Molecule exhibits characteristics consistent with the limonoid class"