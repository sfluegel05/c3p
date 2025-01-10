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

    # Check for triterpene characteristic patterns (complex patterns may vary, example given)
    # Here, use a broad pattern for a polycyclic triterpene
    triterpene_pattern = Chem.MolFromSmarts('C1CCC2(CCC3C(CCC4C3CCC5C4CCCC5)C2)C1')
    if not mol.HasSubstructMatch(triterpene_pattern):
        return False, "No triterpene backbone found"
    
    # Check for glycosidic linkage patterns
    # Example pattern for sugar moiety attachment -O[C@H] (glycosidic ether linkage)
    glycoside_pattern = Chem.MolFromSmarts('O[C@H]')
    if not mol.HasSubstructMatch(glycoside_pattern):
        return False, "No glycosidic linkage found"

    # Check for additional sugar moiety patterns (to confirm presence of polysaccharides)
    sugar_pattern = Chem.MolFromSmarts('C1OC(O)C(O)C(O)C1')
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if len(sugar_matches) < 1:
        return False, "No sugar moieties detected"

    # Verifications for common triterpenoid features can be added as needed here
    return True, "Structure matches triterpenoid saponin characteristics"

# Example execution
# print(is_triterpenoid_saponin("put-a-smiles-string-here"))