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

    # Molecular Weight check - cucurbitacins are triterpenoids
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
      return False, f"Molecular weight {mol_wt} is too low (less than 400 Da)"
    
    # Check for tetracyclic core using a more specific SMARTS pattern
    # This pattern looks for a series of connected rings forming the core
    # It does NOT require them to be fused directly
    # Note: The pattern below identifies key carbon atoms forming the core.
    core_pattern = Chem.MolFromSmarts("[C]1~[C]~[C]~[C]2~[C]~[C]~[C]3~[C]~[C]~[C]4~[C]~[C]~[C]~[C]1~[C]2~[C]34")

    if not mol.HasSubstructMatch(core_pattern):
        return False, "Tetracyclic core not found"
          
    # If all checks passed, we consider it a cucurbitacin
    return True, "Likely a cucurbitacin based on structure features"