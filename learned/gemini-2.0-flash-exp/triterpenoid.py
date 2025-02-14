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
    ring_pattern = Chem.MolFromSmarts("[R]1[R][R][R]1") #minimal ring with 4 atoms
    ring_matches = mol.GetSubstructMatches(ring_pattern)
    if len(ring_matches) < 4:
      return False, f"Too few rings, got {len(ring_matches)}, need at least 4."

    # Check carbon count
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 27 or c_count > 33 :
       return False, f"Carbon count out of range {c_count}. Triterpenoids must have ~30 C."
        
    # Check oxygen count
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 1 or o_count > 10:
        return False, f"Oxygen count out of range {o_count}, requires 1 to 10 oxygens."

    # Check methyl groups
    methyl_pattern = Chem.MolFromSmarts("[CH3]")
    methyl_matches = mol.GetSubstructMatches(methyl_pattern)
    if len(methyl_matches) < 4 :
        return False, f"Too few methyl groups, got {len(methyl_matches)}, needs at least 4."


    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400 or mol_wt > 800:
        return False, f"Molecular weight out of range {mol_wt}, must be between 400 and 800."


    return True, "Meets criteria for a triterpenoid (multiple rings, appropriate C count, O count, methyl count, and molecular weight)"