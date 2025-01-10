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

    # General pattern for proteinogenic amino acids
    pattern = Chem.MolFromSmarts("[NX3,NX4][C@@H]([*!$(*=*)])[CH,N,O,S,P]*C(=O)O")
    
    # Check if the molecule matches the generalized amino acid pattern
    if mol.HasSubstructMatch(pattern):
        return True, "Matches generalized proteinogenic amino acid structure"

    return False, "Does not match the generalized structure of a proteinogenic amino acid"


# Example test for the function
smiles_list = [
    'OC(=O)[C@@H](N)CC1=C(C(=C(C(=C1[2H])[2H])[2H])[2H])[2H]',  # L-phenylalanine-d5
    'N[C@@H](Cc1c[nH]cn1)C(O)=O',  # L-histidine
    'C[C@H](N)C(O)=O',  # L-alanine
    'CC(C)[C@H](N)C(O)=O'  # L-valine
]

for smi in smiles_list:
    is_amino_acid, reason = is_proteinogenic_amino_acid(smi)
    print(f"SMILES: {smi}, Is Proteinogenic Amino Acid: {is_amino_acid}, Reason: {reason}")