"""
Classifies: CHEBI:16366 anthocyanidin cation
"""
"""
Classifies: anthocyanidin cation
"""

from rdkit import Chem

def is_anthocyanidin_cation(smiles: str):
    """
    Determines if a molecule is an anthocyanidin cation based on its SMILES string.
    An anthocyanidin cation is an aglycone of anthocyanin cation; they are oxygenated derivatives 
    of flavylium (2-phenylchromenylium).
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anthocyanidin cation, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for positive charge
    charge = mol.GetFormalCharge()
    if charge <= 0:
        return False, "Molecule is not a cation"

    # Define SMARTS pattern for flavylium core (2-phenylchromenylium cation)
    # Flavylium core: fused tricyclic ring system with a positively charged oxygen atom
    flavylium_pattern = Chem.MolFromSmarts('c1ccccc1c2cc[o+]c3ccccc23')
    if not mol.HasSubstructMatch(flavylium_pattern):
        return False, "Flavylium core not found"

    # Check for oxygenated substituents (hydroxyl or methoxy groups) on aromatic rings
    # Phenolic OH groups and methoxy groups attached to aromatic carbons
    oxygenated = False
    oxy_substituent_pattern = Chem.MolFromSmarts('[cH]O')  # Phenolic OH
    methoxy_pattern = Chem.MolFromSmarts('[cH]OC')  # Methoxy group
    if mol.HasSubstructMatch(oxy_substituent_pattern):
        oxygenated = True
    elif mol.HasSubstructMatch(methoxy_pattern):
        oxygenated = True
    else:
        return False, "No oxygenated substituents found on aromatic rings"

    # Check for absence of sugar moieties (aglycone)
    # Look for glycosidic bonds: C-O-C between carbons not in rings
    sugar_pattern = Chem.MolFromSmarts('[#6;R0]-O-[#6;R0]')
    if mol.HasSubstructMatch(sugar_pattern):
        return False, "Sugar moieties detected (glycosidic bonds present)"

    return True, "Molecule is an anthocyanidin cation"