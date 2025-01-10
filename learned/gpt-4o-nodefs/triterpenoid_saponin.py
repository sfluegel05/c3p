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

    # Generalized triterpenoid-related backbones: Oleanane, Ursolane, etc.
    triterpene_patterns = [
        Chem.MolFromSmarts('C1CCC2(CC[C@]3(CC[C@]4(C)C3[C@H](CC[C@@H]5[C@](C)(CCC2[C@H]1)C4=O)C)C)C'),
        Chem.MolFromSmarts('C1CCC2(CC[C@]3(CC[C@]4(C)[C@]3(C)CCC2[C@H]1)C)C'),
    ]

    if not any(mol.HasSubstructMatch(pattern) for pattern in triterpene_patterns):
        return False, "No triterpene backbone found"

    # More generalized glycoside attachment
    glycoside_pattern = Chem.MolFromSmarts('O[C@H]')
    if not mol.HasSubstructMatch(glycoside_pattern):
        return False, "No glycosidic linkage found"

    # Flexible pattern for sugar moieties (supporting round a wide variety of configurations)
    polysaccharide_pattern = Chem.MolFromSmarts('C1OC(O[C@@H])C(O)C(O)C1')
    polysaccharide_matches = mol.GetSubstructMatches(polysaccharide_pattern)
    if len(polysaccharide_matches) < 1:
        return False, "No polysaccharide sugar moieties detected"

    # Additional checks for atom types and counts could be included here
    return True, "Structure matches triterpenoid saponin characteristics"

# Example execution
# print(is_triterpenoid_saponin("put-a-smiles-string-here"))