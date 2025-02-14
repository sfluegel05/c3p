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

    # Define SMARTS pattern for alkanethiol (-SH group attached to a sp3 carbon that is attached to at least one carbon)
    # Handles cases with double bonds: [SH][C]([#6])[#6] or [SH][CX4]([#6])[#6] which covers the case C(S)C
    #  and [SH][C]=[C]
    alkanethiol_pattern1 = Chem.MolFromSmarts("[SH][C]([#6])[#6]")
    alkanethiol_pattern2 = Chem.MolFromSmarts("[SH][CX4]([#6])")
    alkanethiol_pattern3 = Chem.MolFromSmarts("[SH][C]=[C]")

    # Check if the molecule matches the pattern
    if mol.HasSubstructMatch(alkanethiol_pattern1) or mol.HasSubstructMatch(alkanethiol_pattern2) or mol.HasSubstructMatch(alkanethiol_pattern3):
        
        matches = []
        if mol.HasSubstructMatch(alkanethiol_pattern1):
            matches.extend(mol.GetSubstructMatches(alkanethiol_pattern1))
        if mol.HasSubstructMatch(alkanethiol_pattern2):
            matches.extend(mol.GetSubstructMatches(alkanethiol_pattern2))
        if mol.HasSubstructMatch(alkanethiol_pattern3):
            matches.extend(mol.GetSubstructMatches(alkanethiol_pattern3))
        
        for match in matches:
            sulfur_atom = mol.GetAtomWithIdx(match[0])
            carbon_atom = mol.GetAtomWithIdx(match[1])
            
            if carbon_atom.IsInRing():
                continue;

            return True, "Contains a sulfanyl group (-SH) attached to an alkyl group."
    

    return False, "Does not contain a sulfanyl group (-SH) attached to an alkyl group."