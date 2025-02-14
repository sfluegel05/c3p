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

    # Define the myo-inositol stereochemistry pattern
    # Myo-inositol configuration follows [C@H] or [C@@H] pattern correctly alternating
    myo_inositol_pattern = Chem.MolFromSmarts("C1([C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@@H]1O)")
    if not mol.HasSubstructMatch(myo_inositol_pattern):
        return False, "Molecule does not match the myo-inositol stereochemistry"

    # Check for the presence of phosphate groups (consider both -P(O)(O)=O and -P(O)(=O)-O forms)
    phosphate_pattern1 = Chem.MolFromSmarts("OP(=O)(O)O")
    phosphate_pattern2 = Chem.MolFromSmarts("OP(=O)(O)O")  # alternative phosphate visualization
    phosphate_count1 = len(mol.GetSubstructMatches(phosphate_pattern1))
    phosphate_count2 = len(mol.GetSubstructMatches(phosphate_pattern2))
    if phosphate_count1 + phosphate_count2 == 0:
        return False, "No phosphate groups found"

    return True, "Molecule matches the myo-inositol stereochemistry with phosphate groups attached"