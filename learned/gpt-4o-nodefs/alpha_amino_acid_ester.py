"""
Classifies: CHEBI:46874 alpha-amino acid ester
"""
from rdkit import Chem

def is_alpha_amino_acid_ester(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid ester based on its SMILES string.

    An alpha-amino acid ester typically has:
    - An alpha-carbon connecting to an amino group.
    - A carboxyl group (-COO-) esterified.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is an alpha-amino acid ester, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern: capture alpha-amino acid ester structures
    # Allow both linear and cyclic structures
    # Pattern: alpha-carbon [CH3 or CH] connected to amino [N] and esterified carboxyl [C(=O)O]
    alpha_amino_ester_pattern1 = Chem.MolFromSmarts('[CH1,CH2]([NX3])-[CX3](=O)-[OX2][CX4,CX3]')
    alpha_amino_ester_pattern2 = Chem.MolFromSmarts('[CX4]([NX3])[CX3](=O)-[O][C]')
    alpha_amino_ester_pattern3 = Chem.MolFromSmarts('C(=O)[OX2][CX4,CX3][NX3]')

    if mol.HasSubstructMatch(alpha_amino_ester_pattern1) or mol.HasSubstructMatch(alpha_amino_ester_pattern2) or mol.HasSubstructMatch(alpha_amino_ester_pattern3):
        return True, "Molecule contains an alpha-amino acid ester structure"
    else:
        return False, "No alpha-amino acid ester pattern found"

# Example usage
example_smiles = "COC(=O)CN"
result, reason = is_alpha_amino_acid_ester(example_smiles)
print(f"Classification result: {result}, Reason: {reason}")