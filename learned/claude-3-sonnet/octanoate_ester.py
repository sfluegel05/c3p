"""
Classifies: CHEBI:87657 octanoate ester
"""
"""
Classifies: CHEBI:35997 octanoate ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_octanoate_ester(smiles: str):
    """
    Determines if a molecule is an octanoate ester based on its SMILES string.
    An octanoate ester is any fatty acid ester in which the carboxylic acid component is octanoic acid (caprylic acid).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an octanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define patterns for octanoic acid and ester bonds
    octanoic_acid_pattern = Chem.MolFromSmarts("CCCCCCCC(=O)[OH]")
    ester_pattern = Chem.MolFromSmarts("C(=O)OC")
    
    # Check for octanoic acid and ester substructures
    octanoic_acid_matches = mol.GetSubstructMatches(octanoic_acid_pattern)
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    if not octanoic_acid_matches or not ester_matches:
        return False, "Missing octanoic acid or ester substructure"
    
    # Check if any ester bond is attached to octanoic acid
    for ester_match in ester_matches:
        ester_atom = mol.GetAtomWithIdx(ester_match[1])
        if ester_atom.GetDegree() == 2:  # Exclude false positives like ring systems
            neighbor_atoms = [mol.GetAtomWithIdx(nbr_idx).GetSmarts() for nbr_idx in ester_atom.GetNeighbors()]
            if "CCCCCCCC(=O)[OH]" in neighbor_atoms:
                return True, "Contains octanoic acid ester group"
    
    return False, "No octanoic acid ester group found"