"""
Classifies: CHEBI:26493 quinic acid
"""
from rdkit import Chem

def is_quinic_acid(smiles: str):
    """
    Determines if a molecule is a quinic acid or its derivative based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quinic acid derivative, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string into RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Quinic acid core: cyclohexane with 3 hydroxy groups, 1 carboxylic acid
    quinic_acid_base_pattern = Chem.MolFromSmarts("OC1CC(O)C(O)CC1C(O)=O")
    
    # Check if molecule contains the quinic acid core
    if not mol.HasSubstructMatch(quinic_acid_base_pattern):
        return False, "Mismatch in quinic acid core structure"

    # Check for presence of ester linkages derived from cinnamic acid derivatives (e.g., caffeoyl)
    ester_pattern = Chem.MolFromSmarts("O=C(O)COC=Cc1cc(O)c(O)c1")  # representative ester group pattern
    if mol.HasSubstructMatch(ester_pattern):
        return True, "Matches quinic acid core structure with caffeoyl derivative"

    # Finally, confirm as quinic acid if the core pattern is there and no other complex ester is detected
    return True, "Matches quinic acid core structure"