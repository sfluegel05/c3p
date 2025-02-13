"""
Classifies: CHEBI:15341 beta-D-glucosiduronic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_beta_D_glucosiduronic_acid(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronic acid based on its SMILES string.
    A beta-D-glucosiduronic acid is a glucuronic acid condensed with another substance
    to form a glycosidic bond, with the glucuronic acid in the beta configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-glucosiduronic acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define glucuronic acid SMARTS pattern
    glucuronic_acid_pattern = Chem.MolFromSmarts("[C@H]1(O)[C@@H](O)[C@H](O)[C@@H](O[C@@H]2[C@@H](C(=O)O)O[C@H](CO)[C@@H]2O)[C@H](O)1")

    # Check for glucuronic acid substructure
    matches = mol.GetSubstructMatches(glucuronic_acid_pattern)
    if not matches:
        return False, "No glucuronic acid substructure found"

    # Check for glycosidic bond (O-C-O) between glucuronic acid and aglycone
    for match in matches:
        anom_carbon_idx = match[3]  # Index of anomeric carbon
        anom_carbon = mol.GetAtomWithIdx(anom_carbon_idx)
        if sum(mol.GetAtomWithIdx(n).GetAtomicNum() == 8 for n in anom_carbon.GetNeighbors()) == 2:
            break  # Found glycosidic bond
    else:
        return False, "No glycosidic bond found between glucuronic acid and aglycone"

    # Check for beta configuration
    if mol.GetAtomWithIdx(match[4]).GetChiralTag() != Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
        return False, "Glucuronic acid not in beta configuration"

    return True, "Contains glucuronic acid in beta configuration, connected via glycosidic bond to an aglycone"