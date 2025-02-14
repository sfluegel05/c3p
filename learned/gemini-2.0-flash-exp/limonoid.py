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

    # Check for furan ring with possible substituents
    furan_pattern = Chem.MolFromSmarts("[c1cco[c,C]1]")
    if not mol.HasSubstructMatch(furan_pattern):
        return False, "No furan ring found."

    # Check for a tetracyclic or higher core (at least 4 fused rings)
    # We look for ring systems with at least 4 rings based on connectivity
    num_rings = Chem.GetSSSR(mol)
    if num_rings < 4:
        return False, "No tetracyclic or higher core found."
    
    # Check the overall number of carbons (range)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 25 or carbon_count > 40:
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