"""
Classifies: CHEBI:33856 aromatic amino acid
"""
from rdkit import Chem

def is_aromatic_amino_acid(smiles: str):
    """
    Determines if a molecule is an aromatic amino acid based on its SMILES string.
    An aromatic amino acid is characterized by the presence of a connected amino group, 
    a carboxylic acid group, and an aromatic ring system.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule matches aromatic amino acid characteristics, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for key functional features in aromatic amino acids
    amino_pattern = Chem.MolFromSmarts("[NX3;H2,H1,H0+0]")  # NH2 or NH group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O;H0,H1-]")  # COOH or similar
    aromatic_system_pattern = Chem.MolFromSmarts("a1aaaaa1")  # Phenyl-like full aromatic ring

    # Assert presence of an amino group directly connected to a carbon that forms the backbone (typical for amino acids)
    if not mol.HasSubstructMatch(amino_pattern):
        return False, "No amino group found"
    
    # Assert presence of a carboxylic group
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid or related group found"
    
    # Ensure an aromatic ring system is part of the structure, ensuring full ring presence
    if not mol.HasSubstructMatch(aromatic_system_pattern):
        return False, "No full aromatic ring system found"

    return True, "Contains amino group, carboxylic acid, and aromatic ring system"

# Testing with provided examples
examples = [
    "Nc1ccc(cn1)C(O)=O", # 6-aminonicotinic acid
    "Nc1cc(O)cc(c1)C(O)=O", # 3-amino-5-hydroxybenzoic acid
    "N[C@@H](Cc1c[nH]cn1)C(O)=O", # L-histidine
    # Add more as needed for testing
]

for example in examples:
    print(is_aromatic_amino_acid(example))