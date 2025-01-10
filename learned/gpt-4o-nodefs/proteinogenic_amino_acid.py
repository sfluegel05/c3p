"""
Classifies: CHEBI:83813 proteinogenic amino acid
"""
from rdkit import Chem

def is_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a proteinogenic amino acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a proteinogenic amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Amino acid pattern: basic form with alpha carbon, amino group, and carboxylic acid
    amino_acid_pattern = Chem.MolFromSmarts("N[C@@H](C)C(=O)O")  # Simplified pattern
    if mol.HasSubstructMatch(amino_acid_pattern):
        return True, "Matches basic proteinogenic amino acid structure"

    return False, "Does not match the structure of a proteinogenic amino acid"


# Testing the function with example SMILES strings
smiles_list = [
    'OC(=O)[C@@H](N)CC1=C(C(=C(C(=C1[2H])[2H])[2H])[2H])[2H]',  # L-phenylalanine-d5
    'N[C@@H](Cc1c[nH]cn1)C(O)=O',  # L-histidine
    'C[C@H](N)C(O)=O',  # L-alanine
    'CC(C)[C@H](N)C(O)=O'  # L-valine
]

for smiles in smiles_list:
    is_amino_acid, reason = is_proteinogenic_amino_acid(smiles)
    print(f"SMILES: {smiles}, Is Proteinogenic Amino Acid: {is_amino_acid}, Reason: {reason}")