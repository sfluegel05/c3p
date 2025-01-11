"""
Classifies: CHEBI:134045 polychlorinated dibenzodioxines and related compounds
"""
from rdkit import Chem

def is_polychlorinated_dibenzodioxines_and_related_compounds(smiles: str):
    """
    Determines if a molecule belongs to the class of polychlorinated dibenzodioxines and related compounds.
    
    This includes checking for chlorinated dibenzodioxins, dibenzofurans, and biphenyls structures
    with an emphasis on polychlorination.
    
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

    # Count chlorine and bromine atoms for broader halogenation checking
    halogen_count = sum(atom.GetSymbol() in ['Cl', 'Br'] for atom in mol.GetAtoms())
    if halogen_count < 4:
        return False, "Less than four halogen atoms"

    # Broader SMARTS patterns for variability
    dibenzodioxin_pattern = Chem.MolFromSmarts("c1cc2oc3ccccc3Oc2cc1")
    dibenzofuran_pattern = Chem.MolFromSmarts("c1cc2oc3ccccc3c2cc1")
    biphenyl_pattern = Chem.MolFromSmarts("c1ccccc1-c2ccccc2")

    # Check for substructure matches with broader chlorination requirements
    if mol.HasSubstructMatch(dibenzodioxin_pattern):
        return True, "Polychlorinated dibenzodioxin structure found"
    elif mol.HasSubstructMatch(dibenzofuran_pattern):
        return True, "Polychlorinated dibenzofuran structure found"
    elif mol.HasSubstructMatch(biphenyl_pattern):
        return True, "Polyhalogenated biphenyl structure found"

    return False, "Does not match the specific known structures for this class"