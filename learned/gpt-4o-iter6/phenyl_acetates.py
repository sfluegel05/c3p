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

    # Define a more accurate SMARTS pattern for acetoxy group (OC(=O)C) directly connected to an aromatic ring
    phenyl_acetate_pattern = Chem.MolFromSmarts("[cH]1c([OH])ccc1OC(=O)C")
    
    # Check if the molecule matches the phenotype acetate pattern
    if mol.HasSubstructMatch(phenyl_acetate_pattern):
        return True, "Molecule is a phenyl acetate with correctly positioned acetate ester linked to an aromatic carbon"

    # Extended check for phenyl-like substructures – capture a broader set of similar aromatic ester links
    extended_pattern = Chem.MolFromSmarts("[c]1ccc([OH])c1OC(=O)C")
    if mol.HasSubstructMatch(extended_pattern):
        return True, "Molecule exhibits an extended phenyl acetate structure"

    return False, "Molecule does not exhibit key phenyl acetate structural traits"