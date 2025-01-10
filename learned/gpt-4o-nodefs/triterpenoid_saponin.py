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

    # Define triterpene backbones: Oleanane, Ursolane
    triterpene_patterns = [
        Chem.MolFromSmarts('C1CCC2(CC[C@]3(CC[C@]4(C)C3CCC[C@]4(C2)C1)C)C'),  # Oleanane-like
        Chem.MolFromSmarts('C1CCC2(CC[C@]3(CC[C@]4(C)C3CCC[C@@]4(C2)C1)C)C'),  # Ursane-like
    ]

    if not any(mol.HasSubstructMatch(pattern) for pattern in triterpene_patterns):
        return False, "No triterpene backbone found"

    # Detect glycosidic linkages (general pattern)
    glycoside_pattern = Chem.MolFromSmarts('[CX4:1][OX2:2][C]')
    if not mol.HasSubstructMatch(glycoside_pattern):
        return False, "No glycosidic linkage found"

    # Detect sugar moieties (e.g., glucose, rhamnose); flexible to various ring forms
    sugar_patterns = [
        Chem.MolFromSmarts('C1O[C@H](C)CO[C@H]1O'),  # Glucose-like
        Chem.MolFromSmarts('C1O[C@H](C)[C@@H](O)[C@H](O)[C@H]1O')  # Rhamnose-like
    ]

    sugar_count = 0
    for pattern in sugar_patterns:
        sugar_count += len(mol.GetSubstructMatches(pattern))
        
    if sugar_count < 1:
        return False, "No sugar moieties detected"

    return True, "Structure matches triterpenoid saponin characteristics"

# An example execution is not provided due to complexity, but test with real SMILES strings.