"""
Classifies: CHEBI:35436 D-glucoside
"""
"""
Classifies: CHEBI:XXXXX D-glucoside (any glucoside derived from D-glucose)
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_D_glucoside(smiles: str):
    """
    Determines if a molecule is a D-glucoside based on its SMILES string.
    A D-glucoside contains a D-glucose moiety connected via glycosidic bond.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define D-glucose core pattern (pyranose form with correct stereochemistry)
    # This pattern matches beta-D-glucopyranose with glycosidic oxygen
    glucose_pattern = Chem.MolFromSmarts(
        "[C@@H]1([C@H]([C@@H]([C@H]([C@H](O1)CO)O)O)O)O-O-*"
    )
    
    # Find matches for D-glucose core with glycosidic bond
    matches = mol.GetSubstructMatches(glucose_pattern)
    if not matches:
        # Try alpha-D configuration alternative
        alpha_glucose_pattern = Chem.MolFromSmarts(
            "[C@H]1([C@H]([C@@H]([C@H]([C@H](O1)CO)O)O)O)O-O-*"
        )
        if not mol.GetSubstructMatches(alpha_glucose_pattern):
            return False, "No D-glucose core with glycosidic bond found"

    # Verify the glucose unit has exactly 5 hydroxyl groups
    hydroxyl_count = sum(1 for atom in mol.GetAtoms() 
                        if atom.GetAtomicNum() == 8 and 
                        atom.GetTotalNumHs() >= 1 and 
                        atom.GetDegree() == 2)
    if hydroxyl_count < 5:
        return False, "Insufficient hydroxyl groups for glucose unit"

    # Check molecular formula contains C6H10O5 (glucose core)
    formula = rdMolDescriptors.CalcMolFormula(mol)
    if "C6H10O5" not in formula:
        return False, "Missing glucose core formula components"

    return True, "Contains D-glucose moiety with glycosidic bond"