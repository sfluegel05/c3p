"""
Classifies: CHEBI:87658 decanoate ester
"""
"""
Classifies: CHEBI:75840 decanoate ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_decanoate_ester(smiles: str):
    """
    Determines if a molecule is a decanoate ester based on its SMILES string.
    A decanoate ester is a fatty acid ester resulting from the formal condensation of the carboxy group of decanoic acid (capric acid) with the hydroxy group of an alcohol or phenol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a decanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the ester group (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester group found"

    # Check if any ester group is attached to a 10-carbon chain
    decanoate_pattern = Chem.MolFromSmarts("[CX4H2,CX3H1,CX2H0][CX4H2,CX3H1,CX2H0][CX4H2,CX3H1,CX2H0][CX4H2,CX3H1,CX2H0][CX4H2,CX3H1,CX2H0][CX4H2,CX3H1,CX2H0][CX4H2,CX3H1,CX2H0][CX4H2,CX3H1,CX2H0][CX4H2,CX3H1,CX2H0][CX4H3,CX3H2,CX2H1]")
    for match in ester_matches:
        # Get the atom connected to the ester oxygen
        ester_oxygen = match[0]
        ester_carbon = match[1]
        # Get the atom connected to the ester carbon (should be part of the 10-carbon chain)
        ester_carbon_neighbors = mol.GetAtomWithIdx(ester_carbon).GetNeighbors()
        for neighbor in ester_carbon_neighbors:
            if neighbor.GetIdx() != ester_oxygen:
                # Check if this neighbor is part of a 10-carbon chain
                if mol.HasSubstructMatch(decanoate_pattern):
                    return True, "Contains a decanoate ester group (10-carbon chain with ester linkage)"

    return False, "No decanoate ester group found (10-carbon chain with ester linkage)"