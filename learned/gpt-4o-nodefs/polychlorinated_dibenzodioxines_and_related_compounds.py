"""
Classifies: CHEBI:134045 polychlorinated dibenzodioxines and related compounds
"""
from rdkit import Chem

def is_polychlorinated_dibenzodioxines_and_related_compounds(smiles: str):
    """
    Determines if a molecule belongs to the class of polychlorinated dibenzodioxines and related compounds.
    
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

    # Check for at least four chlorine atoms for targeting polychlorinated classes
    cl_count = sum(atom.GetSymbol() == 'Cl' for atom in mol.GetAtoms())
    if cl_count < 4:
        return False, "Less than four chlorine atoms"
    
    # Check specifically for the dibenzodioxin structure with chlorines
    dibenzodioxin_pattern = Chem.MolFromSmarts("c1cc2oc3ccccc3Oc2cc1")
    dibenzodioxin_match = mol.HasSubstructMatch(dibenzodioxin_pattern)
    
    # Check specifically for the dibenzofuran structure with chlorines
    dibenzofuran_pattern = Chem.MolFromSmarts("c1cc2oc3ccccc3c2cc1")
    dibenzofuran_match = mol.HasSubstructMatch(dibenzofuran_pattern)
    
    # Check specifically for polyhalogenated biphenyl structures
    biphenyl_pattern = Chem.MolFromSmarts("c1ccccc1-c2ccccc2")
    biphenyl_match = mol.HasSubstructMatch(biphenyl_pattern)

    if dibenzodioxin_match and cl_count >= 4:
        return True, "Polychlorinated dibenzodioxin structure found"
    elif dibenzofuran_match and cl_count >= 4:
        return True, "Polychlorinated dibenzofuran structure found"
    elif biphenyl_match and cl_count >= 4:
        return True, "Polyhalogenated biphenyl structure found"

    return False, "Does not match the specific known structures for this class"