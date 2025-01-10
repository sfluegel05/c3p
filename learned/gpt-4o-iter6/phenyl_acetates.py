"""
Classifies: CHEBI:140310 phenyl acetates
"""
from rdkit import Chem

def is_phenyl_acetates(smiles: str):
    """
    Determines if a molecule is a phenyl acetate based on its SMILES string.
    A phenyl acetate is defined as an acetate ester obtained by formal condensation of acetic acid with
    the hydroxy group of any phenol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a phenyl acetate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define more specific pattern for acetoxy group connected to an aromatic ring (phenyl group)
    # The pattern first ensures an aromatic ring is connected to an ester group OC(=O)C
    phenyl_acetate_pattern = Chem.MolFromSmarts("c1ccccc1[O][C](=O)C")
    
    # Check if the molecule matches the phenyl acetate pattern
    if mol.HasSubstructMatch(phenyl_acetate_pattern):
        return True, "Molecule is a phenyl acetate with correctly positioned acetate ester linked to an aromatic carbon"

    return False, "Molecule does not exhibit key phenyl acetate structural traits"