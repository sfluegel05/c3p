from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_biflavonoid(smiles: str):
    """
    Determines if a molecule is a biflavonoid based on the presence of two flavonoid units
    connected by a single atom or bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a biflavonoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for minimum molecular weight (roughly 2x flavonoid MW)
    mw = Chem.Descriptors.ExactMolWt(mol)
    if mw < 450:  # Approximate minimum MW for biflavonoids
        return False, "Molecular weight too low for biflavonoid"

    # Basic flavonoid core patterns
    patterns = [
        # Flavone core
        'c1cccc(c1)-c1cc(=O)c2c(O)cccc2o1',
        # Flavan core
        'c1cccc(c1)C1Oc2cccc(O)c2CC1',
        # Flavanone core
        'c1cccc(c1)C1CC(=O)c2c(O)cccc2O1'
    ]

    total_matches = 0
    for pattern in patterns:
        patt = Chem.MolFromSmarts(pattern)
        if patt is not None:
            matches = mol.GetSubstructMatches(patt)
            total_matches += len(matches)

    # Count aromatic rings
    aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    
    # Count oxygen atoms
    num_oxygens = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8])

    # Check for characteristic ketone/hydroxyl groups
    ketone_pattern = Chem.MolFromSmarts('C(=O)')
    hydroxyl_pattern = Chem.MolFromSmarts('[OH]')
    
    ketone_count = len(mol.GetSubstructMatches(ketone_pattern)) if ketone_pattern else 0
    hydroxyl_count = len(mol.GetSubstructMatches(hydroxyl_pattern)) if hydroxyl_pattern else 0

    # Criteria for classification
    if aromatic_rings < 4:
        return False, "Insufficient number of aromatic rings for biflavonoid"
    
    if num_oxygens < 6:
        return False, "Insufficient number of oxygen atoms for biflavonoid"
        
    if ketone_count + hydroxyl_count < 4:
        return False, "Insufficient number of ketone/hydroxyl groups"

    # Check for presence of flavonoid-like structures
    if total_matches < 1:
        return False, "Does not contain required flavonoid core structures"

    # Additional check for common biflavonoid features
    if (aromatic_rings >= 4 and 
        num_oxygens >= 8 and 
        (ketone_count >= 2 or hydroxyl_count >= 4) and 
        mw >= 450 and
        total_matches >= 1):
        return True, "Contains features consistent with biflavonoid structure"

    return False, "Does not meet structural requirements for biflavonoid classification"
# Pr=0.7285714285714285
# Recall=0.9272727272727272