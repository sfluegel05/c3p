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

    # Check for the presence of at least one furan ring
    furan_pattern = Chem.MolFromSmarts("c1ccoc1")
    if not mol.HasSubstructMatch(furan_pattern):
        return False, "No furan ring identified, which is typical for limonoid"

    # Inspect for a potential limonoid core pattern
    limonoid_core_pattern = Chem.MolFromSmarts("C1CCC2CCC3C4CCC5=C3C2=C1C=5")  # Simplified core pattern
    if not mol.HasSubstructMatch(limonoid_core_pattern):
        return False, "Lacks core structural features typical of limonoids"

    # Evaluate the extent of oxygenation through number of oxygen atoms
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 5:
        return False, "Insufficient oxygen atoms for high oxygenation typical of limonoids"
 
    # Consider molecular weight as a supportive factor
    mol_wt = Descriptors.MolWt(mol)
    # While generally realistic, avoid overstrict limitations
    if mol_wt < 400 or mol_wt > 1200:
        return False, "Molecular weight falls outside typical range but still possible limonoid"

    return True, "Molecule exhibits characteristics consistent with the limonoid class"