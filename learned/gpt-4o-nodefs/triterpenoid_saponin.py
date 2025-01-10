"""
Classifies: CHEBI:61778 triterpenoid saponin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_triterpenoid_saponin(smiles: str):
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define triterpene backbones: Oleanane, Ursane, Dammarane
    triterpene_patterns = [
        Chem.MolFromSmarts('C1CCC2(CC[C@]3(CC[C@]4(C)C3CCC[C@]4(C2)C1)C)C'),  # Oleanane-like
        Chem.MolFromSmarts('C1CCC2(CC[C@]3(CC[C@]4(C)C3CCC[C@]4(C2)C1)C)C'),  # Ursane-like
        Chem.MolFromSmarts('C1CCC2(C3CCC4(C3CCC34CCCC1)C)C4CC2'),  # Dammarane-like
    ]

    if not any(mol.HasSubstructMatch(pattern) for pattern in triterpene_patterns):
        return False, "No triterpenoid backbone found"

    # Detect glycosidic linkages (particular focus on ether/ester bonds)
    glycoside_pattern = Chem.MolFromSmarts('O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1')  # Example for glucose linkage
    if not mol.HasSubstructMatch(glycoside_pattern):
        return False, "No glycosidic linkage found"

    # Detect common sugar moieties (e.g., glucose, rhamnose, xylose)
    sugar_patterns = [
        Chem.MolFromSmarts('C1O[C@H](CO)[C@H](O)[C@@H](O)[C@H]1O'),  # Glucose
        Chem.MolFromSmarts('C1O[C@H](C)[C@@H](O)[C@H](O)[C@H]1O'),  # Rhamnose
        Chem.MolFromSmarts('C1O[C@H](C)[C@@H](O)COC1=O'),          # Xylose (linear)
    ]

    sugar_count = 0
    for pattern in sugar_patterns:
        sugar_count += len(mol.GetSubstructMatches(pattern))
        
    if sugar_count < 1:
        return False, "No sugar moieties detected"

    return True, "Structure matches triterpenoid saponin characteristics"

# This function can be tested by invoking it with valid SMILES strings of known triterpenoid saponins.