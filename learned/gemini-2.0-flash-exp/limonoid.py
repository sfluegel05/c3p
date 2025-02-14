"""
Classifies: CHEBI:39434 limonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_limonoid(smiles: str):
    """
    Determines if a molecule is a limonoid based on its SMILES string.
    Limonoids are triterpenoids with a characteristic 4,4,8-trimethyl-17-furanylsteroid skeleton
    and high degree of oxygenation.

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

    # Check for furan ring
    furan_pattern = Chem.MolFromSmarts("c1ccoc1")
    if not mol.HasSubstructMatch(furan_pattern):
        return False, "No furan ring found."

     # Check for a triterpenoid core (4 fused rings) with some flexibility.
    # This pattern is not perfect, but it captures the basic tetracyclic skeleton
    core_pattern = Chem.MolFromSmarts("[C]1[C]([C])[C]2[C]([C])[C]([C])[C]3[C]([C])([C])([C])[C]([C])([C])[C]([C])([C])([C])C([C])([C])C3C2C1")
    if not mol.HasSubstructMatch(core_pattern):
         return False, "No tetracyclic triterpenoid core found."


    # Check the overall number of carbons (range)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 20 or carbon_count > 40:
        return False, f"Number of carbons is outside triterpenoid range, found {carbon_count}"

    # Check for high oxygenation (at least 6 oxygens)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 6:
        return False, f"Too few oxygen atoms, found {oxygen_count}, need at least 6"
    
    # Check Molecular Weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 350:
        return False, f"Molecular weight too low for a limonoid, {mol_wt}"


    # If all criteria pass return true
    return True, "Matches limonoid criteria: triterpenoid with furan, tetracyclic core, within carbon range, and high oxygenation."