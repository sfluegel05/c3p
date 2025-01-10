"""
Classifies: CHEBI:48039 dihydroflavonols
"""
"""
Classifies: CHEBI:XXXXX dihydroflavonols
Dihydroflavonols are hydroxyflavanones with a hydroxy group at position 3 of the heterocyclic ring.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_dihydroflavonols(smiles: str):
    """
    Determines if a molecule is a dihydroflavonol based on its SMILES string.
    A dihydroflavonol is a hydroxyflavanone with a hydroxy group at position 3 of the heterocyclic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dihydroflavonol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more flexible flavanone backbone pattern
    flavanone_pattern = Chem.MolFromSmarts("[O]=C1C(=O)C=C2C(=C1)OC=C2")
    if not mol.HasSubstructMatch(flavanone_pattern):
        return False, "No flavanone backbone found"

    # Define the pattern for a hydroxyl group at position 3 of the heterocyclic ring
    hydroxyl_at_3_pattern = Chem.MolFromSmarts("[OH]C1C(=O)C=C2C(=C1)OC=C2")
    if not mol.HasSubstructMatch(hydroxyl_at_3_pattern):
        return False, "No hydroxyl group at position 3 of the heterocyclic ring"

    # Check for additional hydroxyl groups (optional, but common in dihydroflavonols)
    hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() == 1)
    if hydroxyl_count < 2:
        return False, "Insufficient hydroxyl groups for a typical dihydroflavonol"

    # Check molecular weight (dihydroflavonols typically have MW > 200 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, "Molecular weight too low for a dihydroflavonol"

    return True, "Contains flavanone backbone with a hydroxyl group at position 3 of the heterocyclic ring"