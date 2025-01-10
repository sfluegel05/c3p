"""
Classifies: CHEBI:132155 hydroxynaphthoquinone
"""
"""
Classifies: CHEBI:XXXXX hydroxynaphthoquinone
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_hydroxynaphthoquinone(smiles: str):
    """
    Determines if a molecule is a hydroxynaphthoquinone based on its SMILES string.
    A hydroxynaphthoquinone is a naphthoquinone with at least one hydroxy group attached.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxynaphthoquinone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define naphthoquinone core patterns (1,4-naphthoquinone and 1,2-naphthoquinone)
    naphthoquinone_pattern_1 = Chem.MolFromSmarts("c1ccc2C(=O)C=CC(=O)c2c1")
    naphthoquinone_pattern_2 = Chem.MolFromSmarts("c1ccc2C(=O)C(=O)C=Cc2c1")
    
    # Check for naphthoquinone core
    if not (mol.HasSubstructMatch(naphthoquinone_pattern_1) or mol.HasSubstructMatch(naphthoquinone_pattern_2)):
        return False, "No naphthoquinone core found"

    # Define the hydroxy group pattern
    hydroxy_pattern = Chem.MolFromSmarts("[OH]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    if len(hydroxy_matches) == 0:
        return False, "No hydroxy group found"

    # Get the atoms in the naphthoquinone core
    naphthoquinone_atoms = set()
    if mol.HasSubstructMatch(naphthoquinone_pattern_1):
        naphthoquinone_atoms.update(mol.GetSubstructMatch(naphthoquinone_pattern_1))
    if mol.HasSubstructMatch(naphthoquinone_pattern_2):
        naphthoquinone_atoms.update(mol.GetSubstructMatch(naphthoquinone_pattern_2))

    # Check if at least one hydroxy group is directly attached to the naphthoquinone core or on a side chain directly connected to the core
    for match in hydroxy_matches:
        hydroxy_atom = match[0]
        for neighbor in mol.GetAtomWithIdx(hydroxy_atom).GetNeighbors():
            if neighbor.GetIdx() in naphthoquinone_atoms:
                return True, "Contains a naphthoquinone core with at least one hydroxy group attached"
            # Check if the hydroxy group is on a side chain directly connected to the core
            for second_neighbor in neighbor.GetNeighbors():
                if second_neighbor.GetIdx() in naphthoquinone_atoms:
                    return True, "Contains a naphthoquinone core with at least one hydroxy group on a side chain directly connected to the core"

    return False, "Hydroxy group not directly attached to the naphthoquinone core or on a side chain directly connected to the core"