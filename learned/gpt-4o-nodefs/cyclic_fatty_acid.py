"""
Classifies: CHEBI:59238 cyclic fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdmolops

def is_cyclic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a cyclic fatty acid based on its SMILES string.
    Cyclic fatty acids generally contain a cycle within a long aliphatic chain
    and a carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyclic fatty acid, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for the presence of at least one ring
    if rdmolops.GetSSSR(mol) == 0:
        return False, "No cyclic structure detected"

    # Check for sufficient length of carbon chain
    carbon_chain_length = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_chain_length < 8:
        return False, f"Carbon chain too short: {carbon_chain_length} carbons"

    # Broaden the recognition of cyclic substructures
    substructure_patterns = [
        Chem.MolFromSmarts("c1ccoc1"),       # Furan ring
        Chem.MolFromSmarts("O=C1OC=CC1"),    # Lactone ring
        Chem.MolFromSmarts("C1OC1"),         # Epoxide
        Chem.MolFromSmarts("C1CCC1"),        # Cyclobutane
        Chem.MolFromSmarts("C1CCOC1"),       # Tetrahydrofuran
        Chem.MolFromSmarts("C1CC=CC1"),      # Cyclohexane variant
        Chem.MolFromSmarts("O=C1CCC(=O)O1"), # Larger Lactone variant
    ]
    
    for pattern in substructure_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, f"Contains a ring structure ({Chem.MolToSmarts(pattern)}) typical of cyclic fatty acids"

    return False, "Does not contain a recognized cyclic fatty acid substructure"