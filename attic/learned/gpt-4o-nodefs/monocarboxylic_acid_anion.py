"""
Classifies: CHEBI:35757 monocarboxylic acid anion
"""
from rdkit import Chem

def is_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a monocarboxylic acid anion based on its SMILES string.
    A monocarboxylic acid anion typically has one deprotonated carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monocarboxylic acid anion, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for exactly one carboxylate group [O-]C(=O)
    carboxylate_pattern = Chem.MolFromSmarts("[O-]C(=O)")
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    
    if len(carboxylate_matches) == 1:
        return True, "Contains exactly one carboxylate group"
    else:
        return False, f"Contains {len(carboxylate_matches)} carboxylate groups, expected exactly 1"

# Example SMILES that belong to this class
test_smiles = [
    "[H][C@]12C[C@@H](O)[C@]([H])(\\C=C\\[C@@H](O)C(C)CC#CC)[C@@]1([H])C1=C(O2)C(CCCC([O-])=O)=CC=C1",
    "C(CCCCCC/C=C/[C@@H](C[C@H]1[C@@H](CCCCC)O1)O)([O-])=O",
    "[H]C(=OC=C[O-])c1ccccc1C([O-])=O"
]

for smiles in test_smiles:
    result, reason = is_monocarboxylic_acid_anion(smiles)
    print(f"SMILES: {smiles}\nResult: {result}, Reason: {reason}\n")