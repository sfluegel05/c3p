"""
Classifies: CHEBI:15841 polypeptide
"""
"""
Classifies: CHEBI:16670 polypeptide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polypeptide(smiles: str):
    """
    Determines if a molecule is a polypeptide (10+ amino acid residues) based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polypeptide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Calculate the number of amide bonds
    num_amide = rdMolDescriptors.CalcNumAmideBonds(mol)
    
    # For linear peptides: residues = amide bonds + 1
    # For cyclic peptides: residues = amide bonds
    # Check if either case gives >=10 residues
    if num_amide >= 9:  # Covers linear (9+1=10) or cyclic (>=10)
        return True, f"Contains at least {num_amide} amide bonds, indicating 10 or more residues"
    else:
        return False, f"Only {num_amide} amide bonds, insufficient for 10 residues"