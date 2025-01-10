"""
Classifies: CHEBI:134045 polychlorinated dibenzodioxines and related compounds
"""
from rdkit import Chem

def is_polychlorinated_dibenzodioxines_and_related_compounds(smiles: str):
    """
    Determines if a molecule belongs to the class of polychlorinated dibenzodioxines and related compounds.
    
    This includes checking for chlorinated dibenzodioxins, dibenzofurans, and biphenyls structures
    with an emphasis on polychlorination, while excluding other functional group variations.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule fits the class description, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count halogen atoms
    halogen_count = sum(atom.GetSymbol() in ['Cl', 'Br'] for atom in mol.GetAtoms())
    
    # Classification relevant thresholds
    if halogen_count < 4:
        return False, "Less than four halogen atoms"
    
    # SMARTS patterns for structural bases
    dibenzodioxin_pattern = Chem.MolFromSmarts("c1cc2oc3ccccc3Oc2cc1")
    dibenzofuran_pattern = Chem.MolFromSmarts("c1cc2oc3ccccc3c2cc1")
    biphenyl_pattern = Chem.MolFromSmarts("c1ccccc1-c2ccccc2")
    
    # Check if molecule is substructure with required polychlorination without extra functional groups
    if mol.HasSubstructMatch(dibenzodioxin_pattern):
        if all(neigh.GetSymbol() in ['C', 'Cl', 'Br'] for atom in mol.GetAtoms() for neigh in atom.GetNeighbors()):
            return True, "Polychlorinated dibenzodioxin structure found"
    elif mol.HasSubstructMatch(dibenzofuran_pattern):
        if all(neigh.GetSymbol() in ['C', 'Cl', 'Br'] for atom in mol.GetAtoms() for neigh in atom.GetNeighbors()):
            return True, "Polychlorinated dibenzofuran structure found"
    elif mol.HasSubstructMatch(biphenyl_pattern):
        if all(neigh.GetSymbol() in ['C', 'Cl', 'Br'] for atom in mol.GetAtoms() for neigh in atom.GetNeighbors()):
            return True, "Polyhalogenated biphenyl structure found"

    return False, "Does not match the specific known structures for this class"