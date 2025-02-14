"""
Classifies: CHEBI:26493 quinic acid
"""
"""
Classifies: CHEBI:18198 quinic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_quinic_acid(smiles: str):
    """
    Determines if a molecule is a quinic acid based on its SMILES string.
    A quinic acid is a cyclitol carboxylic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quinic acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Quinic acid core pattern (cyclohexane ring with carboxylic acid and multiple hydroxyls)
    quinic_acid_pattern = Chem.MolFromSmarts("[C@H]1([C@H]([C@@H]([C@@H]([C@H](C1)C(=O)O)O)O)O)O")

    # Common substituent patterns
    substituent_pattern = Chem.MolFromSmarts("[O;X2]-[C;X3]=[O;X1]")

    # Check for quinic acid core
    if not mol.HasSubstructMatch(quinic_acid_pattern):
        return False, "No quinic acid core found"

    # Check for substituents
    sub_matches = mol.GetSubstructMatches(substituent_pattern)
    num_subs = len(sub_matches)

    # Check molecular descriptors
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    # Classification criteria
    if num_subs == 0 and mol_wt > 192 and mol_wt < 196 and c_count == 7 and o_count == 6:
        return True, "Unsubstituted quinic acid"
    elif num_subs > 0 and mol_wt > 300 and c_count > 10 and o_count > 6:
        return True, f"Quinic acid with {num_subs} substituents"

    return False, "Not a quinic acid derivative"