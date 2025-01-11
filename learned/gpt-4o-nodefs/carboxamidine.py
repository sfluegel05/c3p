"""
Classifies: CHEBI:35359 carboxamidine
"""
"""
Classifies: carboxamidine
"""
from rdkit import Chem

def is_carboxamidine(smiles: str):
    """
    Determines if a molecule contains a carboxamidine group based on its SMILES string.
    A carboxamidine is broadly defined as R-N=C(NR)NR or similar structures where R can vary widely.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains a carboxamidine group, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Enhanced array of potential SMARTS patterns for carboxamidine groups
    carboxamidine_patterns = [
        Chem.MolFromSmarts("N=C(N)N"),            # Basic carboxamidine form
        Chem.MolFromSmarts("[NX3][CX3]=[NX2][NX3]"),  # Varied substitutions on nitrogen and carbon
        Chem.MolFromSmarts("N=C(N)[#6]"),         # Carboxamidine with broad carbon substituents
        Chem.MolFromSmarts("N=C(N)N-C=O"),        # Carboxamidine adjacent to carbonyl, as seen in complex entities
        Chem.MolFromSmarts("[CX3]=[NX2]-[NX3]"),  # Broader match to include aromatic nitrogen atoms
        Chem.MolFromSmarts("N=C(N)(N)[#6]"),      # Adding compound feature for nitrogen-bound functionalities
        Chem.MolFromSmarts("N=C(N)N-[#6]"),       # Presence of further carbon based moieties through N
    ]

    # Check each pattern
    for pattern in carboxamidine_patterns:
        if pattern is None:
            continue
        if mol.HasSubstructMatch(pattern):
            return True, "Contains carboxamidine group matching the structure"

    return False, "Does not contain recognizable carboxamidine groups"