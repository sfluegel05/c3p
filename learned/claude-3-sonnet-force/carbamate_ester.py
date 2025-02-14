"""
Classifies: CHEBI:23003 carbamate ester
"""
"""
Classifies: CHEBI:35669 carbamate ester
'Any ester of carbamic acid or its N-substituted derivatives.'
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_carbamate_ester(smiles: str):
    """
    Determines if a molecule is a carbamate ester based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbamate ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carbamic acid or N-substituted carbamic acid substructure
    carbamic_acid_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])[NX3]")
    carbamic_acid_matches = mol.GetSubstructMatches(carbamic_acid_pattern)

    if not carbamic_acid_matches:
        return False, "No carbamic acid or N-substituted carbamic acid substructure found"

    # Check for ester linkage on the oxygen atom
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    if not ester_matches:
        return False, "Carbamic acid oxygen not esterified"

    # Allow various substituents on the nitrogen atom
    allowed_substituents = ["C", "c", "O", "N", "S", "P", "F", "Cl", "Br", "I"]
    n_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]

    for n_atom in n_atoms:
        neighbors = [mol.GetAtomWithIdx(idx) for idx in n_atom.GetNeighbors()]
        substituents = [neighbor.GetSymbol() for neighbor in neighbors]
        if any(sub not in allowed_substituents for sub in substituents):
            return False, "Disallowed substituent on carbamate nitrogen"

    return True, "Contains carbamate ester substructure with allowed substituents"