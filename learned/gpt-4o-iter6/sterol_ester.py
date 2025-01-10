"""
Classifies: CHEBI:35915 sterol ester
"""
from rdkit import Chem

def is_sterol_ester(smiles: str):
    """
    Determines if a molecule is a sterol ester based on its SMILES string.
    A sterol ester is obtained by formal condensation of 
    the carboxy group of any carboxylic acid with the 3-hydroxy group of a sterol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a general sterol backbone pattern (fused four-ring core)
    sterol_pattern = Chem.MolFromSmarts("C1CCC2C3CCC4CCCC(C)(C3)O2(C3=CC=CC=C3)C1")
    if not mol.HasSubstructMatch(sterol_pattern):
        return False, "No sterol backbone found"

    # Define ester linkage pattern (ester linked to 3-hydroxy group of sterol)
    # "[#6](=O)[O][#6]" defines a generic ester linkage
    # "3-hydroxy" is assumed in sterols, so check connectivity
    ester_linkage_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.HasSubstructMatch(ester_linkage_pattern)
    if not ester_matches:
        return False, "No ester linkage found"

    # Further verify that the ester linkage is connected to the sterol backbone
    for match in mol.GetSubstructMatches(sterol_pattern):
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() == 'C':
                neighbors = atom.GetNeighbors()
                for neighbor in neighbors:
                    if mol.HasSubstructMatch(Chem.MolFromSmarts("[OH]")) and neighbor.GetSymbol() == 'O':
                        return True, "Contains sterol backbone with appropriate ester linkage"

    return False, "No connection between sterol backbone and ester linkage found"