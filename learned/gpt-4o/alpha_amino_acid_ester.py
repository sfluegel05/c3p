"""
Classifies: CHEBI:46874 alpha-amino acid ester
"""
from rdkit import Chem

def is_alpha_amino_acid_ester(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid ester based on its SMILES string.
    An alpha-amino acid ester arises from the condensation of an alpha-amino acid with an alcohol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino acid ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Pattern for recognizing alpha-amino acid esters
    # This SMARTS matches: backbone of an alpha-amino acid ester, NH-CH(R)-C(O)-O-R'
    # Enhanced consideration for alpha-carbon, ester linkage and flexible attachment
    pattern = Chem.MolFromSmarts("[NX3]-[CH](-[CX4,CX3,CX2])-[CX3](=O)-O-[C]")
    
    if pattern is None:
        return None, "Pattern not defined correctly"

    if mol.HasSubstructMatch(pattern):
        return True, "Contains the alpha-amino acid ester functional group"

    return False, "Does not contain the alpha-amino acid ester functional group"