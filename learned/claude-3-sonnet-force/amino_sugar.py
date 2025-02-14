"""
Classifies: CHEBI:28963 amino sugar
"""
"""
Classifies: CHEBI:26651 amino sugar

An amino sugar is any sugar having one or more alcoholic hydroxy groups replaced by
substituted or unsubstituted amino groups.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_amino_sugar(smiles: str):
    """
    Determines if a molecule is an amino sugar based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amino sugar, False otherwise
        str: Reason for classification
    """
    try:
        # Parse SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "Invalid SMILES string"

        # Check for sugar backbone
        sugar_backbone_pattern = Chem.MolFromSmarts("[OX2]C[CX4][CX4][CX4][CX4][CX4][OX2]")
        if not mol.HasSubstructMatch(sugar_backbone_pattern):
            return False, "No sugar backbone found"

        # Check for amino groups replacing hydroxy groups
        amino_pattern = Chem.MolFromSmarts("[NX3][CX4][OX2]")
        amino_matches = mol.GetSubstructMatches(amino_pattern)
        if not amino_matches:
            return False, "No amino groups replacing hydroxy groups"

        # Check for remaining hydroxy groups
        hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
        hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
        if not hydroxy_matches:
            return False, "No remaining hydroxy groups found"

        # Negative check: exclude sugar alcohols
        sugar_alcohol_pattern = Chem.MolFromSmarts("[CX4][OX2H][CX4][OX2H]")
        if mol.HasSubstructMatch(sugar_alcohol_pattern):
            return False, "Molecule appears to be a sugar alcohol, not an amino sugar"

        return True, "Contains a sugar backbone with one or more amino groups replacing hydroxy groups"

    except Exception as e:
        return False, f"An error occurred: {str(e)}"