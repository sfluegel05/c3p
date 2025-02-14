"""
Classifies: CHEBI:32952 amine
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

def is_amine(smiles: str):
    """
    Determines if a molecule is an amine based on its SMILES string.
    An amine is derived from ammonia by replacing one or more hydrogen atoms with hydrocarbyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an amine, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a query for a nitrogen atom with up to three single-bonded carbons (amine N)
    amine_query = rdqueries.AtomNumEqualsQueryAtom(7)  # Nitrogen atom
    amine_query.SetDegree(1, 3)  # Allow between 1 and 3 neighbors

    # Exclude nitrogen atoms in specific functional groups
    non_amine_patterns = [
        Chem.MolFromSmarts("N-C(=O)"),  # Amide
        Chem.MolFromSmarts("N(=O)"),    # Nitro
    ]
    
    is_any_amine = False
    for atom in mol.GetAtoms():
        if amine_query.MatchAtom(atom):
            # Check if this nitrogen is part of any excluded group
            if any(mol.HasSubstructMatch(pattern) for pattern in non_amine_patterns if mol.HasSubstructMatch(pattern, recursive=False)):
                continue
            # If not part of excluded groups and it's an amine pattern
            is_any_amine = True
            break
    
    if is_any_amine:
        return True, "Molecule contains an amine group"
    return False, "No amine group detected"