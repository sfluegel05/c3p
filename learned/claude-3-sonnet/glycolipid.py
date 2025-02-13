"""
Classifies: CHEBI:33563 glycolipid
"""
"""
Classifies: CHEBI:18194 glycolipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glycolipid(smiles: str):
    """
    Determines if a molecule is a glycolipid based on its SMILES string.
    A glycolipid has a carbohydrate part linked to a lipid part.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycolipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carbohydrate part - detect common carbohydrate rings
    glyco_smarts = ['OC1OC(CO)C(O)C1O', 'OC1OC(O)C(O)C(O)C1O'] # glucose, galactose
    glyco_found = False
    for smarts in glyco_smarts:
        glyco_pattern = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(glyco_pattern):
            glyco_found = True
            break
    if not glyco_found:
        return False, "No carbohydrate part found"

    # Look for lipid part - long aliphatic chains
    lipid_smarts = '[C;H3][C;H2][C;H2][C;H2][C;H2][C;H2][C;H2][C;H2]' # at least 8 carbons
    lipid_pattern = Chem.MolFromSmarts(lipid_smarts)
    lipid_matches = mol.GetSubstructMatches(lipid_pattern)
    if not lipid_matches:
        return False, "No lipid part found"

    # Look for glycosidic linkage between carbohydrate and lipid parts
    glyco_linkage_smarts = '[OX2]C[OX2]'
    glyco_linkage_pattern = Chem.MolFromSmarts(glyco_linkage_smarts)
    glyco_linkage_matches = mol.GetSubstructMatches(glyco_linkage_pattern)
    if not glyco_linkage_matches:
        return False, "No glycosidic linkage found between carbohydrate and lipid parts"

    return True, "Contains carbohydrate and lipid parts linked via a glycosidic bond"