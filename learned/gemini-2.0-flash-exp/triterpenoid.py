"""
Classifies: CHEBI:36615 triterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_triterpenoid(smiles: str):
    """
    Determines if a molecule is a triterpenoid based on its SMILES string.
    Triterpenoids are characterized by a 30-carbon skeleton, often with multiple fused rings.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
         bool: True if molecule is a triterpenoid, False otherwise
         str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for a minimum number of rings
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    if num_rings < 2:
      return False, f"Too few rings, got {num_rings}, need at least 2."

    # Check carbon count
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20 or c_count > 80:
       return False, f"Carbon count out of range {c_count}. Triterpenoids must have ~30 C (but can have modifications)."
        
    # Check oxygen count
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 0 or o_count > 25:
        return False, f"Oxygen count out of range {o_count}, requires between 0 to 25 oxygens."

    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 1500:
        return False, f"Molecular weight out of range {mol_wt}, must be between 300 and 1500."


    # Define a SMARTS pattern for the core triterpenoid structure (multiple fused rings with 4 or 5 member rings)
    # This is a very generic pattern and may be too broad, but it's a starting point, and can be made more precise later if needed.
    # I am specifically looking for a 5-membered carbon ring fused to at least two other carbon rings.
    triterpenoid_core_smarts = "[C;R5]1([C;R])([C;R])[C;R]2([C;R])([C;R])[C;R]3([C;R])([C;R])123"
    core_pattern = Chem.MolFromSmarts(triterpenoid_core_smarts)
    if core_pattern is None:
         return False, "Invalid SMARTS pattern."
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Does not match the triterpenoid core structure"


    return True, "Meets criteria for a triterpenoid (multiple rings, appropriate C count, O count, molecular weight, and core structure)."