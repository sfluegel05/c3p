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

    # Define SMARTS pattern for S-C bond (any carbon)
    sc_pattern = Chem.MolFromSmarts("[SH][C]")

    # Define SMARTS pattern for alkyl carbon. This considers all sp3 hybridized carbons.
    alkyl_carbon_pattern = Chem.MolFromSmarts("[CX4]")


    # Check if the molecule has a S-C bond
    if not mol.HasSubstructMatch(sc_pattern):
        return False, "Does not contain a sulfanyl group (-SH) attached to a carbon."
    
    # Get all matches of S-C
    sc_matches = mol.GetSubstructMatches(sc_pattern)
    
    # Verify that at least one carbon bound to S is part of alkyl group
    for match in sc_matches:
        sulfur_atom = mol.GetAtomWithIdx(match[0])
        carbon_atom = mol.GetAtomWithIdx(match[1])

        if carbon_atom.HasSubstructMatch(alkyl_carbon_pattern): # carbon is sp3 hybridized
             
             is_alkyl_chain = False
             for neighbor in carbon_atom.GetNeighbors(): #check that this carbon is attached to other carbons or hydrogens
                 if neighbor.GetAtomicNum() == 6 or neighbor.GetAtomicNum() == 1:
                     is_alkyl_chain = True
                     break
             if is_alkyl_chain:
                return True, "Contains a sulfanyl group (-SH) attached to an alkyl group."
    
    return False, "Does not contain a sulfanyl group (-SH) attached to an alkyl group."