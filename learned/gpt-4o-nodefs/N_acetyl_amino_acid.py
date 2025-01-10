"""
Classifies: CHEBI:21575 N-acetyl-amino acid
"""
from rdkit import Chem

def is_N_acetyl_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acetyl-amino acid based on its SMILES string.
    An N-acetyl-amino acid has an acetyl group attached to the nitrogen of a characteristic amino acid backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple: (bool, str) - True if molecule is an N-acetyl-amino acid, False otherwise along with a reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for acetyl group (CH3-C(=O)) attached to nitrogen
    acetyl_pattern = Chem.MolFromSmarts('CC(=O)N')
    if not mol.HasSubstructMatch(acetyl_pattern):
        return False, "No acetyl group attached to nitrogen found"

    # Verify presence of amino acid carboxylic group
    carboxylate_pattern = Chem.MolFromSmarts('C(=O)[O-]')  # Accounting for the charged conditions too
    carboxylic_pattern = Chem.MolFromSmarts('C(=O)O')
    if not mol.HasSubstructMatch(carboxylic_pattern) and not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "No carboxylic acid group found"

    # Check for the amino acid pattern that includes the backbone being bound directly to an acetylated nitrogen
    # It should look like: N-CH-C(=O)O connecting the acetyl part to the rest
    backbone_pattern = Chem.MolFromSmarts('[NX3]([C;H])[CH1]-C(=O)[O]')
    if not mol.HasSubstructMatch(backbone_pattern):
        return False, "No proper amino acid backbone connected to acetyl group"

    return True, "Contains acetyl group correctly bound to nitrogen of an amino acid"