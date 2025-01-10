"""
Classifies: CHEBI:26255 prenylquinone
"""
"""
Classifies: Prenylquinone
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_prenylquinone(smiles: str):
    """
    Determines if a molecule is a prenylquinone based on its SMILES string.
    A prenylquinone is a quinone substituted by a polyprenyl-derived side-chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prenylquinone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Expand quinone pattern list
    quinone_patterns = [
        Chem.MolFromSmarts("O=C1C=CC(=O)C=C1"),            # Benzoquinone
        Chem.MolFromSmarts("O=C1C=CC(=O)C2=CC(C=CC2)=C1"), # Naphthoquinone
        Chem.MolFromSmarts("O=C1C=C(O)C(=O)C2=CC=CC=C12"), # Anthraquinone
        Chem.MolFromSmarts("O=C1CC(C(=O)C=C1)=O"),         # Hydroxyqu...

    ]

    # Look for any quinone backbone match
    has_quinone = any(mol.HasSubstructMatch(pattern) for pattern in quinone_patterns)
    if not has_quinone:
        return False, "No quinone backbone found"

    # Identify prenyl side-chains; complex patterning for longer chains
    isoprene_unit = Chem.MolFromSmarts("C(=C)CC")
    prenyl_matches = mol.GetSubstructMatches(isoprene_unit)

    # Utilize a heuristic for identifying significantly long prenyl chains (e.g., more than 14 heavy atoms)
    side_chains = [len(match) for match in prenyl_matches]
    if sum(side_chains) < 14:  # Example threshold for length
        return False, "Prenyl side-chain too short"

    return True, "Contains quinone backbone with adequate prenyl side-chain and length"