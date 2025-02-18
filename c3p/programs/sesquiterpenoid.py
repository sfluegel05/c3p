"""
Classifies: CHEBI:26658 sesquiterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a sesquiterpenoid based on its SMILES string.
    Sesquiterpenoids are characterized by a 15-carbon skeleton derived from three isoprene units,
    and can have modifications such as removal of methyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sesquiterpenoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Basic carbon count - must have between 10 and 18 carbons (allows for some modifications)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (10 <= c_count <= 18):
       return False, f"Incorrect number of carbons: {c_count}, must be between 10 and 18."

    # 2. Check for the core carbon skeleton - a 15 carbon chain (acyclic) or a 15 carbon ringed system

    core_pattern_acyclic = Chem.MolFromSmarts("[C;R0]([C;R0])([C;R0])([C;R0])[C;R0][C;R0]([C;R0])([C;R0])([C;R0])[C;R0][C;R0]([C;R0])([C;R0])[C;R0]") # acyclic, 15 carbons
    core_pattern_cyclic = Chem.MolFromSmarts("[C;R]1([C;R])[C;R]2([C;R])[C;R]3([C;R])([C;R]1)[C;R]([C;R]2)[C;R]3") # generalized ringed
    
    found_core = mol.HasSubstructMatch(core_pattern_acyclic) or mol.HasSubstructMatch(core_pattern_cyclic)
    
    if not found_core:
        return False, "No common sesquiterpenoid skeleton found"


    return True, "Contains a sesquiterpenoid core with 10-18 carbons."