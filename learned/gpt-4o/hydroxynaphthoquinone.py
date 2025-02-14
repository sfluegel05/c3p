"""
Classifies: CHEBI:132155 hydroxynaphthoquinone
"""
from rdkit import Chem

def is_hydroxynaphthoquinone(smiles: str):
    """
    Determines if a molecule is a hydroxynaphthoquinone based on its SMILES string.
    A hydroxynaphthoquinone is defined as any naphthoquinone in which the naphthoquinone moiety
    is substituted by at least one hydroxy group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a hydroxynaphthoquinone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # A more flexible SMARTS pattern for naphthoquinone with possible aromatic substitutions
    naphthoquinone_pattern = Chem.MolFromSmarts("C1=CC=C2C(=O)C=CC(=O)C2=C1")
    naphthoquinone_matches = mol.GetSubstructMatches(naphthoquinone_pattern)
    if len(naphthoquinone_matches) < 1:
        return False, "No naphthoquinone core structure found"
    
    # Check for hydroxy groups within the molecule
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    if len(hydroxy_matches) < 1:
        return False, "No hydroxy group found"

    # Ensure hydroxy is part of the same aromatic system as naphthoquinone
    for match in naphthoquinone_matches:
        for hydroxy in hydroxy_matches:
            # If an oxygen atom in the hydroxy group is connected to an aromatic carbon in the core
            if any(mol.GetBondBetweenAtoms(hydroxy[0], atom_idx).IsAromatic() for atom_idx in match):
                return True, "Contains naphthoquinone core with at least one hydroxy group attached"

    return False, "Hydroxy groups are not properly connected to the naphthoquinone core"

# Example Usage:
# smiles_example = "O=C1C=CC(=O)c2cc(O)ccc12"  # Example of a hydroxynaphthoquinone
# result, reason = is_hydroxynaphthoquinone(smiles_example)
# print(result, reason)