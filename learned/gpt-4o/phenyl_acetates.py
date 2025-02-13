"""
Classifies: CHEBI:140310 phenyl acetates
"""
from rdkit import Chem

def is_phenyl_acetates(smiles: str):
    """
    Determines if a molecule is a phenyl acetate based on its SMILES string.
    A phenyl acetate is an acetate ester obtained by formal condensation
    of the carboxy group of acetic acid with the hydroxy group of any phenol.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a phenyl acetate, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the ester linkage pattern (acetate group)
    ester_pattern = Chem.MolFromSmarts("OC(=O)C")

    # Define aromatic ring with attachment site
    aromatic_oxy_pattern = Chem.MolFromSmarts("aOc")

    # Check for ester group in the molecule
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No acetate ester linkage found"
    
    # Check for an aromatic connection with ester oxygen
    if not mol.HasSubstructMatch(aromatic_oxy_pattern):
        return False, "No aromatic structure connected through ester to oxy group"

    # Ensure that both patterns are found in the molecule
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    aro_matches = mol.GetSubstructMatches(aromatic_oxy_pattern)
    
    if ester_matches and aro_matches:
        return True, "Contains a phenyl acetate ester linkage"
    
    return False, "Ester linkage and aromatic structure are missing proper connection"