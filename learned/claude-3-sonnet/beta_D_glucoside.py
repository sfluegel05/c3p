"""
Classifies: CHEBI:22798 beta-D-glucoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_D_glucoside(smiles: str):
    """
    Determines if a molecule is a beta-D-glucoside based on its SMILES string.
    A beta-D-glucoside is a D-glucoside with the anomeric center in the beta configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-glucoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Sanitize molecule and assign stereochemistry
    AllChem.AssignStereochemistry(mol, cleanIt=True, force=True)

    # Look for beta-D-glucose substructure pattern
    beta_d_glucose_pattern = Chem.MolFromSmarts('[C@@H]1([C@H]([C@@H]([C@H]([C@H](O1)CO)O)O)O)O')
    beta_d_glucose_matches = mol.GetSubstructMatches(beta_d_glucose_pattern)
    if not beta_d_glucose_matches:
        return False, "No beta-D-glucose substructure found"

    # Look for glucosidic bond pattern (O-[C@H]1[C@@H]...)
    glucosidic_bond_pattern = Chem.MolFromSmarts('[O-][C@H]1[C@@H]([C@H]([C@@H]([C@H](O1)CO)O)O)O')
    glucosidic_bond_matches = mol.GetSubstructMatches(glucosidic_bond_pattern)
    if not glucosidic_bond_matches:
        return False, "No glucosidic bond found"

    # Count oxygen atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count != 6:
        return False, "Must have exactly 6 oxygens"

    return True, "Contains a beta-D-glucose moiety with a glucosidic bond"