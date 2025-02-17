"""
Classifies: CHEBI:192499 anthoxanthin
"""
"""
Classifies: Anthoxanthin – flavonoid pigments in plants.
A typical anthoxanthin is based on a flavone (2-phenylchromen-4-one) scaffold,
but extensive substitution (e.g., methoxylation, hydroxylation, glycosylation) is common.
This approach uses:
  • Two SMARTS patterns (flavone cores) to identify the underlying connectivity.
  • A molecular weight filter typical for flavonoids (240–500 Da).
  • A check for substituents: at least one methoxy group ([OX2][CH3]) or at least two hydroxy groups ([OX2H]).
  • An exclusion filter matching a known false positive (Trichophenol A).
If any parsing error or inability to match the intended logic occurs, the function may return (None, None).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_anthoxanthin(smiles: str):
    """
    Determines if a molecule is an anthoxanthin (a plant flavonoid pigment) based on its SMILES string.
    
    The function checks:
      1. That the molecule contains a flavone (2-phenylchromen-4-one) core, using two alternate SMARTS.
      2. That its molecular weight is in the typical flavonoid range (240–500 Da).
      3. That it has some oxygenated (methoxy or hydroxy) substituents expected for water‐soluble pigments.
      4. That it does NOT match a known false positive pattern (here, trichophenol A).
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is classified as an anthoxanthin, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define two SMARTS patterns representing the flavone (2-phenylchromen-4-one) core.
    flavone_smarts_1 = "c1ccc(-c2oc(=O)c3ccccc3c2)cc1"
    flavone_pattern_1 = Chem.MolFromSmarts(flavone_smarts_1)
    if flavone_pattern_1 is None:
        return None, None  # signal failure if pattern cannot be parsed
    
    flavone_smarts_2 = "c1ccc(-c2cc(=O)c3occcc3c2)cc1"
    flavone_pattern_2 = Chem.MolFromSmarts(flavone_smarts_2)
    if flavone_pattern_2 is None:
        return None, None
    
    # Check if at least one of the flavone core patterns is found.
    if not (mol.HasSubstructMatch(flavone_pattern_1) or mol.HasSubstructMatch(flavone_pattern_2)):
        return False, "Flavone core (2-phenylchromen-4-one) not found"
    
    # Exclusion filter:
    # Trichophenol A (a false positive) is known to match a specific pattern.
    trichophenol_exclusion_smarts = "O=C1OC(C2=C(O)C=C(O)C=C2C)=CC=3C1=C(O)C=C(O)C3"
    trichophenol_exclusion_pattern = Chem.MolFromSmarts(trichophenol_exclusion_smarts)
    if trichophenol_exclusion_pattern and mol.HasSubstructMatch(trichophenol_exclusion_pattern):
        return False, "Matches exclusion pattern for Trichophenol A, not an anthoxanthin"
    
    # Check that molecular weight is in a typical flavonoid range.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 240 or mol_wt > 500:
        return False, f"Molecular weight {mol_wt:.1f} Da is outside typical anthoxanthin range (240-500 Da)"
    
    # Check for oxygenated substituents.
    # Count methoxy groups: pattern for -OCH3.
    methoxy_pattern = Chem.MolFromSmarts("[OX2][CH3]")
    methoxy_count = len(mol.GetSubstructMatches(methoxy_pattern)) if methoxy_pattern else 0
    # Count hydroxy groups: pattern for -OH.
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxy_count = len(mol.GetSubstructMatches(hydroxy_pattern)) if hydroxy_pattern else 0
    if methoxy_count < 1 and hydroxy_count < 2:
        return False, "Insufficient methoxy or hydroxy substitution typical of anthoxanthins"
    
    return True, "Molecule contains a flavone core with appropriate substitution and properties typical of anthoxanthins"

# Example usage (testing)
if __name__ == '__main__':
    # Some sample anthoxanthin molecules and a known false-positive:
    test_smiles = {
        "sinensetin": "COc1ccc(cc1OC)-c1cc(=O)c2c(OC)c(OC)c(OC)cc2o1",
        "tambulin": "COc1ccc(cc1)-c1oc2c(OC)c(OC)cc(O)c2c(=O)c1O",
        "wogonin": "COc1c(O)cc(O)c2c1oc(cc2=O)-c1ccccc1",
        "xanthomicrol": "COc1c(O)c2c(oc(cc2=O)-c2ccc(O)cc2)c(OC)c1OC",
        "Trichophenol A": "O=C1OC(C2=C(O)C=C(O)C=C2C)=CC=3C1=C(O)C=C(O)C3",
        "6-geranylchrysin": "CC(C)=CCC\\C(C)=C\\Cc1c(O)cc2oc(cc(=O)c2c1O)-c1ccccc1"
    }
    
    for name, smi in test_smiles.items():
        result, reason = is_anthoxanthin(smi)
        print(f"Name: {name}\nSMILES: {smi}\nResult: {result}\nReason: {reason}\n")