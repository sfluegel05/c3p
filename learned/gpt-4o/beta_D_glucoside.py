"""
Classifies: CHEBI:22798 beta-D-glucoside
"""
from rdkit import Chem

def is_beta_D_glucoside(smiles: str):
    """
    Determines if a molecule is a beta-D-glucoside based on its SMILES string.
    A beta-D-glucoside is characterized by a beta-D-glucose moiety linked through an ether bond.
    The stereochemistry must reflect the beta-configuration at the anomeric site.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a beta-D-glucoside, False otherwise.
        str: Reason for classification.
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Construct a comprehensive pattern for beta-D-glucoside
    # - Identify the pyranose ring of glucose
    # - Ensure beta stereochemistry at the anomeric carbon
    pattern_glucoside = Chem.MolFromSmarts("""
        [C@H]1(O[C@H](CO[*])[*])O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O
        """)

    if mol.HasSubstructMatch(pattern_glucoside):
        return True, "Contains beta-D-glucoside substructure with beta-configuration at anomeric center"
    
    return False, "No beta-D-glucoside substructure found"