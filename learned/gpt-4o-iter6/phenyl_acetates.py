"""
Classifies: CHEBI:140310 phenyl acetates
"""
from rdkit import Chem

def is_phenyl_acetates(smiles: str):
    """
    Determines if a molecule is a phenyl acetate based on its SMILES string.
    A phenyl acetate is defined as an acetate ester obtained by formal condensation
    of the carboxy group of acetic acid with the hydroxy group of any phenolic compound.

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

    # Define a SMARTS pattern for acetoxy group (OC(=O)C) directly connected to an aromatic ring
    phenyl_acetate_pattern = Chem.MolFromSmarts("c1ccccc1OC(=O)C")

    # Check if the molecule matches the phenyl acetate pattern
    if mol.HasSubstructMatch(phenyl_acetate_pattern):
        return True, "Molecule is a phenyl acetate with correctly positioned acetate ester linked to an aromatic ring"

    # In this case, the common motifs should be captured in a broader sense
    extended_phenyl_pattern = Chem.MolFromSmarts("c1ccc(cc1)OC(=O)C")
    if mol.HasSubstructMatch(extended_phenyl_pattern):
        return True, "Molecule exhibits an extended phenyl acetate structure"

    return False, "Molecule does not exhibit key phenyl acetate structural traits"