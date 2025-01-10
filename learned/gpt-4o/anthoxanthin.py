"""
Classifies: CHEBI:192499 anthoxanthin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_anthoxanthin(smiles: str):
    """
    Determines if a molecule is likely an anthoxanthin based on its SMILES string.
    Anthoxanthins are flavonoid pigments found in plants characterized by flavonoid cores
    and hydroxyl/methoxy/glycosyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely an anthoxanthin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for flavone backbone pattern
    # Pattern: C=O group conjugated with a 6-membered aromatic ring (flavonoid basic structure)
    flavonoid_pattern = Chem.MolFromSmarts("O=Cc1ccccc1")
    
    if not mol.HasSubstructMatch(flavonoid_pattern):
        return False, "No flavone backbone found"

    # Check for presence of multiple hydroxyl (and/or methoxy) groups
    n_hydroxy = sum(atom.GetSymbol() == 'O' and atom.GetDegree() == 1 for atom in mol.GetAtoms())
    
    if n_hydroxy < 2:
        return False, f"Insufficient number of hydroxy groups, found {n_hydroxy}"
    
    # Check for complexity indicating potential glycosides or other substituents (high rotatable bond count)
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 3:
        return False, "Molecule too simple, likely not an anthoxanthin"

    # Additional checks can include searching for methoxy groups and common glycosidic units
    methoxy_pattern = Chem.MolFromSmarts("CO")
    count_methoxy = len(mol.GetSubstructMatches(methoxy_pattern))
    if count_methoxy == 0:
        return True, "Flavone backbone with hydroxy groups, likely anthoxanthin"

    # Final decision based on traits across known anthoxanthins
    return True, "Contains flavone backbone and characteristic groups of anthoxanthins"