"""
Classifies: CHEBI:25448 myo-inositol phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_myo_inositol_phosphate(smiles: str):
    """
    Determines if a molecule is a myo-inositol phosphate based on its SMILES string.
    A myo-inositol phosphate is an inositol structure with myo-configuration that has phosphate groups attached.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a myo-inositol phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the correct core of myo-inositol (cyclohexane with alternating stereochemistry for hydroxyl groups)
    myo_inositol_core = Chem.MolFromSmarts("C1[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O")
    if not mol.HasSubstructMatch(myo_inositol_core):
        return False, "Molecule does not match the myo-inositol core configuration"

    # Search for general phosphate groups
    phosphate_patterns = [
        Chem.MolFromSmarts("OP(O)(=O)[O-]"),
        Chem.MolFromSmarts("OP(O)(=O)O")
    ]

    phosphate_matches = sum(mol.HasSubstructMatch(phosphate) for phosphate in phosphate_patterns)
    if not phosphate_matches:
        return False, "No phosphate groups found"

    # Ensure there are no long carbon chains which indicate non-inositol components (e.g., fats)
    longest_chain = max(len(chain) for chain in Chem.rdmolops.GetMolFrags(mol, asMols=False))
    if longest_chain > 6:
        return False, "Long carbon chain detected, possibly non-inositol component"

    return True, "Molecule contains myo-inositol core with phosphate groups attached"