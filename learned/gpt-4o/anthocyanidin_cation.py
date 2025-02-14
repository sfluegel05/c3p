"""
Classifies: CHEBI:16366 anthocyanidin cation
"""
from rdkit import Chem

def is_anthocyanidin_cation(smiles: str):
    """
    Determines if a molecule is an anthocyanidin cation based on its SMILES string.
    An anthocyanidin cation is characterized by a flavylium ion core with additional oxygenation.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anthocyanidin cation, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (None, "Invalid SMILES string")

    # SMARTS pattern revised for the flavylium core: a 2-phenylchromenylium structure
    # The pattern targets the chromenylium base with a positive charge on oxygen
    flavylium_pattern = Chem.MolFromSmarts("Oc1cc(O)c2[o+]c(ccc2c1)-c3cc(O)c(O)c3")  
    if not mol.HasSubstructMatch(flavylium_pattern):
        return (False, "No flavylium core found or it is incorrectly structured")
    
    # Assess oxygenation via hydroxyl, ether, or methoxy groups directly
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    methoxy_pattern = Chem.MolFromSmarts("[OCH3]")
    ether_pattern = Chem.MolFromSmarts("[O]([#6])[#6]")

    # Count oxygen-containing groups that typically contribute to substitution
    oxy_group_matches = (mol.GetSubstructMatches(hydroxyl_pattern) +
                         mol.GetSubstructMatches(methoxy_pattern) +
                         mol.GetSubstructMatches(ether_pattern))
    
    # Confirm presence of sufficient oxygenation without specifying exact numbers due to variance
    if len(oxy_group_matches) < 2:
        return (False, "Insufficient oxygenation; expected substitution")

    # Verify positive charge presence, ensuring cationic nature
    positive_charge = any(atom.GetFormalCharge() > 0 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if not positive_charge:
        return (False, "Molecule lacks a positive charge indicating it is not a cation")

    return (True, "Contains a correctly structured flavylium ion core with positive charge and sufficient oxygenation")