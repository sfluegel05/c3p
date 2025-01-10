"""
Classifies: CHEBI:28868 fatty acid anion
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid_anion(smiles: str):
    """
    Determines if a molecule is a fatty acid anion based on its SMILES string.
    A fatty acid anion is characterized by a deprotonated carboxylic acid group (-C([O-])=O)
    and a hydrocarbon chain, which may include double bonds, rings, or functional groups,
    and should have sufficient chain length.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acid anion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Require a carboxylate group
    carboxylate_pattern = Chem.MolFromSmarts("C(=O)[O-]")
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "No carboxylate group (-C([O-])=O) found"

    # Check carbon chain length considering C and O atoms as flexible members of the chain
    carbon_pattern = Chem.MolFromSmarts("[C,c]")
    oxy_pattern = Chem.MolFromSmarts("[O,o]")
    carbons = mol.GetSubstructMatches(carbon_pattern)
    oxygens = mol.GetSubstructMatches(oxy_pattern)

    # Convert to a collection that maintains atom indices and excludes duplicates
    carbon_set = set(carbons) | set(oxygens)

    if len(carbon_set) < 8:  # Increased minimum length flexibility
        return False, "Insufficient linear carbon chain length for fatty acid anion"

    # Allow for rings, but limit complexity
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    if n_rings > 3:  # Allows some rings but limits overall complexity
        return False, "Presence of too many rings, which is atypical for fatty acid anions"

    # Checking allowed functional groups
    allowable_functional_groups = ["OH", "=O", "C=C"]  # basic example patterns typically found in fatty acids
    for fg in allowable_functional_groups:
        fg_pattern = Chem.MolFromSmarts(fg)
        if not mol.HasSubstructMatch(fg_pattern):
            continue

    return True, "Contains a deprotonated carboxylate group with a sufficiently long and possible non-linear carbon chain including allowable functional groups"