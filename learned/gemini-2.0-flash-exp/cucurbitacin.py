"""
Classifies: CHEBI:16219 cucurbitacin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_cucurbitacin(smiles: str):
    """
    Determines if a molecule is a cucurbitacin based on its SMILES string.
    Cucurbitacins are tetracyclic triterpenoids with specific structural characteristics.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cucurbitacin, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count Rings (must be 4 or more)
    ring_info = mol.GetRingInfo()
    n_rings = ring_info.NumRings()
    if n_rings < 4:
        return False, "Less than 4 rings"
    
    # Molecular Weight check - cucurbitacins are triterpenoids
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400 or mol_wt > 1000:
      return False, f"Molecular weight {mol_wt} is out of the typical range (400-1000 Da)"
    
    # Check for tetracyclic core using a more specific SMARTS pattern
    # This pattern looks for a fused 4-ring system
    # Note that the rings must be fused!
    core_pattern = Chem.MolFromSmarts("[C]1:[C]:[C]2:[C]:[C]3:[C]:[C]4:[C]1:[C]2:[C]34") #This requires a fused system
    
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Tetracyclic core not found"
          
    # If all checks passed, we consider it a cucurbitacin
    return True, "Likely a cucurbitacin based on structure features"