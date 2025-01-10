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
    containing a beta-lactam ring fused to a five-membered ring, with various substitutions at positions 3, 4, and 6.

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
    # This pattern matches the core bicyclic structure with a beta-lactam ring
    # and allows for various substitutions at positions 3, 4, and 6
    carbapenem_pattern = Chem.MolFromSmarts("[C@H]12[C@H]([C,c])C(=O)N1C=C[C@H]2[C,c]")
    if not mol.HasSubstructMatch(carbapenem_pattern):
        # Try alternative stereochemistry
        carbapenem_pattern = Chem.MolFromSmarts("[C@@H]12[C@@H]([C,c])C(=O)N1C=C[C@@H]2[C,c]")
        if not mol.HasSubstructMatch(carbapenem_pattern):
            return False, "No carbapenem skeleton found"

    # Check for beta-lactam ring
    beta_lactam_pattern = Chem.MolFromSmarts("C1(=O)NCC1")
    if not mol.HasSubstructMatch(beta_lactam_pattern):
        return False, "No beta-lactam ring found"

    # Check molecular weight range (expanded to include more carbapenems)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 150 or mol_wt > 1000:
        return False, "Molecular weight out of typical carbapenem range"

    # Check for common substituents at positions 3, 4, and 6
    # Position 3: often has a sulfur-containing group
    # Position 4: often has a hydroxyl or substituted carbon
    # Position 6: often has a carboxylate or substituted carbon
    substituent_pattern = Chem.MolFromSmarts("[C@H]12[C@H]([C,c])C(=O)N1C=C[C@H]2[C,c]([O,S])")
    if not mol.HasSubstructMatch(substituent_pattern):
        return False, "Missing characteristic substituents"

    return True, "Contains carbapenem skeleton with characteristic substituents"