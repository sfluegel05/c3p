"""
Classifies: CHEBI:71548 dihydroagarofuran sesquiterpenoid
"""
"""
Classifies: CHEBI:51841 dihydroagarofuran sesquiterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_dihydroagarofuran_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a dihydroagarofuran sesquiterpenoid based on its SMILES string.
    A dihydroagarofuran sesquiterpenoid is a sesquiterpenoid containing a dihydroagarofuran skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dihydroagarofuran sesquiterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for dihydroagarofuran core
    dihydroagarofuran_core = Chem.MolFromSmarts("[C@H]1[C@@]2([C@H](O)[C@@H]([C@H](O)[C@@H]1O)[C@@H]2O)O")
    if not mol.HasSubstructMatch(dihydroagarofuran_core):
        return False, "No dihydroagarofuran core found"

    # Check for sesquiterpenoid backbone
    sesquiterpenoid_backbone = Chem.MolFromSmarts("[C@H]1[C@H]2[C@@H]([C@@H]([C@H]([C@@H]1O)O)O)[C@@H]2O")
    if not mol.HasSubstructMatch(sesquiterpenoid_backbone):
        return False, "Sesquiterpenoid backbone not found"

    # Check for typical substituents and modifications
    has_acetyl_groups = any(atom.GetSmarts() == "[#6]-[#6]-[#8]-[#6]" for atom in mol.GetAtoms())
    has_benzoyl_groups = any(atom.GetSmarts() == "[#6]-[#6]-[#8]-[#6]-[#6]" for atom in mol.GetAtoms())

    # Sesquiterpenoid with dihydroagarofuran core and typical substituents/modifications
    reason = "Contains the dihydroagarofuran core and sesquiterpenoid backbone"
    if has_acetyl_groups:
        reason += ", with acetyl groups"
    if has_benzoyl_groups:
        reason += ", with benzoyl groups"

    return True, reason