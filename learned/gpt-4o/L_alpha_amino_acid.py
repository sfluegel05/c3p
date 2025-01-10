"""
Classifies: CHEBI:15705 L-alpha-amino acid
"""
from rdkit import Chem

def is_L_alpha_amino_acid(smiles: str) -> (bool, str):
    """
    Determines if a molecule is an L-alpha-amino acid based on its SMILES string. 
    L-alpha-amino acids have an alpha carbon with L-configuration attached to an amino group, 
    a carboxyl group, and a side chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an L-alpha-amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for L-alpha-amino acid pattern: C attached to N, C, and a carboxylic group
    l_alpha_amino_pattern = Chem.MolFromSmarts("[C@@H](N)(C(O)=O)")
    if not mol.HasSubstructMatch(l_alpha_amino_pattern):
        return False, "Does not match the L-alpha-amino acid pattern"
        
    return True, "Matches the L-alpha-amino acid pattern with L configuration at alpha carbon"

# Examples for testing
smiles_list = [
    "[H][C@@]1(N\C([C@@H](CCC(O)=O)[C@@H]1=)=CC1=CC(=O)c2c(C)c(Cc3[nH]c(CC4NC(=O)C(C)=C4C=C)c(C)c3CC)[nH]c12)C(O)=O",
    "C1(=CNC2=C1C=CC=C2)C([C@@H](C(=O)O)N)O",  # 3-hydroxy-L-tryptophan
    "O=C(O)[C@@H](N)CCCCCSC",  # L-trihomomethionine
    # Add other example SMILES strings here...
]

for smiles in smiles_list:
    result, reason = is_L_alpha_amino_acid(smiles)
    print(f"SMILES: {smiles}\nResult: {result}, Reason: {reason}\n")