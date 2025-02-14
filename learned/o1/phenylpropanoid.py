"""
Classifies: CHEBI:26004 phenylpropanoid
"""
"""
Classifies: Phenylpropanoid
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phenylpropanoid(smiles: str):
    """
    Determines if a molecule is a phenylpropanoid based on its SMILES string.
    A phenylpropanoid is any organic aromatic compound with a structure based on a phenylpropane skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phenylpropanoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for phenylpropane backbone: benzene ring connected to a three-carbon chain
    phenylpropane_pattern = Chem.MolFromSmarts("c1ccccc1CCC")
    if mol.HasSubstructMatch(phenylpropane_pattern):
        return True, "Contains phenylpropane skeleton"

    # Extended patterns to account for derivatives (e.g., double bonds, hydroxyl groups)
    # Pattern for cinnamyl alcohol and similar structures
    cinnamyl_pattern = Chem.MolFromSmarts("c1ccccc1C=CC")
    if mol.HasSubstructMatch(cinnamyl_pattern):
        return True, "Contains cinnamyl-like structure"

    # Check for flavonoid backbone (phenylchromanone)
    flavonoid_pattern = Chem.MolFromSmarts("c1cc(-c2ccc3oc(=O)ccc3c2)ccc1")
    if mol.HasSubstructMatch(flavonoid_pattern):
        return True, "Contains flavonoid backbone"

    # Check for coumarin structure
    coumarin_pattern = Chem.MolFromSmarts("O=C1C=CC2=CC=CC=C2O1")
    if mol.HasSubstructMatch(coumarin_pattern):
        return True, "Contains coumarin structure"

    # Check for lignan structure (dimers of phenylpropanoids)
    lignan_pattern = Chem.MolFromSmarts("c1cc(Cc2ccc(O)cc2)ccc1")
    if mol.HasSubstructMatch(lignan_pattern):
        return True, "Contains lignan structure"

    # Check for stilbenoid structure (e.g., resveratrol)
    stilbenoid_pattern = Chem.MolFromSmarts("c1ccc(cc1)/C=C/c2ccccc2")
    if mol.HasSubstructMatch(stilbenoid_pattern):
        return True, "Contains stilbenoid structure"

    # Check for phenylpropenoid units with modifications
    phenylpropenoid_pattern = Chem.MolFromSmarts("c1ccccc1C=CC(=O)O")
    if mol.HasSubstructMatch(phenylpropenoid_pattern):
        return True, "Contains phenylpropenoid unit with modifications"

    # Check for precursors of lignin (coniferyl alcohol derivatives)
    lignin_precursor_pattern = Chem.MolFromSmarts("c1cc(O)c(C=CCO)cc1")
    if mol.HasSubstructMatch(lignin_precursor_pattern):
        return True, "Contains lignin precursor structure"

    # Check for other small phenolic molecules based on phenylpropane
    phenolic_pattern = Chem.MolFromSmarts("c1cc(O)ccc1CCC")
    if mol.HasSubstructMatch(phenolic_pattern):
        return True, "Contains phenolic phenylpropane structure"

    # If none of the patterns matched
    return False, "Does not contain phenylpropanoid skeleton"