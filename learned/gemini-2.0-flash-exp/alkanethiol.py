"""
Classifies: CHEBI:47908 alkanethiol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_alkanethiol(smiles: str):
    """
    Determines if a molecule is an alkanethiol based on its SMILES string.
    An alkanethiol is a compound in which a sulfanyl group, -SH, is attached to an alkyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkanethiol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the core SMARTS pattern for alkanethiol (-SH group attached to an alkyl C atom)
    alkanethiol_pattern = Chem.MolFromSmarts("[SH][CX4]")

    # Check if the molecule matches the pattern
    if mol.HasSubstructMatch(alkanethiol_pattern):
       
        # Check if one the carbons attached to the S is part of an alkane group
        match = mol.GetSubstructMatches(alkanethiol_pattern)
        for item in match:
            sulfur_atom = mol.GetAtomWithIdx(item[0])
            carbon_atom = mol.GetAtomWithIdx(item[1])
           
            
        return True, "Contains a sulfanyl group (-SH) attached to an alkyl group."

    return False, "Does not contain a sulfanyl group (-SH) attached to an alkyl group."