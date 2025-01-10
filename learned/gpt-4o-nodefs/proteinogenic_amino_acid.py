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
    
    # Parse the SMILES string into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a pattern for a generic proteinogenic amino acid
    # Central carbon (Cα) with appropriate groups
    amino_acid_pattern = Chem.MolFromSmarts(
        "[NX3;H2,H1;!$(NC=O)]-[C@@;!$([C;D2]=O)]-[CX3](=[OX1])[O;H1,0]"  # Amino group, chiral center @ Cα, and carboxylic acid
    )
    
    # Run a substructure search with this pattern
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "Does not match core proteinogenic amino acid structure"
    
    # If found matches to the structural pattern that often characterize the proteinogenic amino acids
    return True, "Matches refined structure expected of proteinogenic amino acids"

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