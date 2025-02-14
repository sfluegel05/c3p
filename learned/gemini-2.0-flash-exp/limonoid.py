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

    # Check for furan ring attached to the core structure (more specific)
    # The [c] is important to avoid matching furan rings outside the core
    furan_pattern = Chem.MolFromSmarts("[c]1[c][o][c][c]1[C]")
    if furan_pattern is None or not mol.HasSubstructMatch(furan_pattern):
        return False, "No furan ring attached to core found."

     # Check for a 4,4,8-trimethyl pattern
    trimethyl_pattern = Chem.MolFromSmarts("C[C](C)([C])")
    if trimethyl_pattern is None or not mol.HasSubstructMatch(trimethyl_pattern):
         return False, "No 4,4,8-trimethyl pattern found."
         
    # Check for tetracyclic or higher core
    num_rings = Chem.GetSSSR(mol)
    if num_rings < 4:
       return False, "Less than 4 rings found"
    
    # Check the overall number of carbons (range)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 25 or carbon_count > 40:
         return False, f"Carbon count outside triterpenoid range: {carbon_count}"


    # Check for high oxygenation (at least 4 oxygens)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 4:
        return False, f"Too few oxygen atoms, found {oxygen_count}, need at least 4"

    # Molecular Weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, f"Molecular weight too low: {mol_wt}"

    # If all criteria pass return true
    return True, "Matches limonoid criteria: triterpenoid with furan, tetracyclic core, within carbon range, and high oxygenation."