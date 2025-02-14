"""
Classifies: CHEBI:23849 diterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_diterpenoid(smiles: str):
    """
    Determines if a molecule is a diterpenoid based on its SMILES string.
    Diterpenoids have a C20 skeleton (or modified) and are built from isoprene units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diterpenoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
     # Check carbon count (allow for some variation around 20 carbons)
    if c_count < 15 or c_count > 40:  # Adjusted range to be more lenient on modifications.
         return False, f"Carbon count is {c_count}, not in the typical diterpenoid range (15-40)"


    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 500 :
         return False, f"Molecular weight of {mol_wt} is not in the typical diterpenoid range (250-500)"
    


    # Check for ring systems
    n_rings = rdMolDescriptors.CalcNumRings(mol)

    if n_rings < 2: #Most diterpenoids should have at least 2 rings.
         return False, f"Diterpenoids should have at least 2 rings, this one has {n_rings}"
    
    # More specific isoprene unit pattern - looking for connected branches
    # This is a more specific SMARTS pattern for two isoprene units connected together
    isoprene_pattern = Chem.MolFromSmarts("[CX4](C)([CX4])[CX4]~[CX4]~[CX4](C)([CX4])[CX4]") # Two connected branched isoprene units

    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)

    if len(isoprene_matches) < 1:  # Require at least one pair of connected isoprenes
        return False, "Too few connected isoprene units detected (less than 1)"
    
    

    return True, "Likely a diterpenoid based on carbon count, ring system, molecular weight and isoprene units"