"""
Classifies: CHEBI:46633 carbapenems
"""
"""
Classifies: CHEBI:60853 carbapenem
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_carbapenems(smiles: str):
    """
    Determines if a molecule is a carbapenem based on its SMILES string.
    A carbapenem is a beta-lactam antibiotic with a carbapenem skeleton, which is a bicyclic structure
    containing a beta-lactam ring fused to a five-membered ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbapenem, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more general carbapenem skeleton pattern
    # This pattern matches the bicyclic structure with a beta-lactam ring
    carbapenem_pattern = Chem.MolFromSmarts("[C@H]12[C@H](C)C(=O)N1C=C[C@H]2C")
    if not mol.HasSubstructMatch(carbapenem_pattern):
        # Try alternative stereochemistry
        carbapenem_pattern = Chem.MolFromSmarts("[C@@H]12[C@@H](C)C(=O)N1C=C[C@@H]2C")
        if not mol.HasSubstructMatch(carbapenem_pattern):
            return False, "No carbapenem skeleton found"

    # Check for beta-lactam ring (more general pattern)
    beta_lactam_pattern = Chem.MolFromSmarts("[C@H]1[C@H](C)C(=O)N1")
    if not mol.HasSubstructMatch(beta_lactam_pattern):
        # Try alternative stereochemistry
        beta_lactam_pattern = Chem.MolFromSmarts("[C@@H]1[C@@H](C)C(=O)N1")
        if not mol.HasSubstructMatch(beta_lactam_pattern):
            return False, "No beta-lactam ring found"

    # Check for common features (not strict positions)
    # At least one of these should be present:
    # - Sulfur-containing group
    # - Carboxyl group
    # - Hydroxyl group
    sulfur_pattern = Chem.MolFromSmarts("[SX2]")
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    
    has_sulfur = mol.HasSubstructMatch(sulfur_pattern)
    has_carboxyl = mol.HasSubstructMatch(carboxyl_pattern)
    has_hydroxyl = mol.HasSubstructMatch(hydroxyl_pattern)
    
    if not (has_sulfur or has_carboxyl or has_hydroxyl):
        return False, "Missing common carbapenem features (sulfur, carboxyl, or hydroxyl groups)"

    # Check molecular weight range (typical for carbapenems)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 600:
        return False, "Molecular weight out of typical carbapenem range"

    return True, "Contains carbapenem skeleton with common features"