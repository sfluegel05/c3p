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
    An anthocyanidin cation is an aglycon of anthocyanin cation; they are oxygenated derivatives 
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
    
    # Check for positively charged oxygen atom in the flavylium core
    # Define the flavylium core pattern based on flavylium cation SMILES
    flavylium_core = Chem.MolFromSmarts('c1ccc(cc1)C2=CC=[O+]C3=CC=CC=C23')
    if not mol.HasSubstructMatch(flavylium_core):
        return False, "Flavylium core not found"

    # Check for oxygenated substituents (e.g., hydroxyl or methoxy groups) on aromatic rings
    # Phenolic OH groups and methoxy groups attached to aromatic carbons
    oxy_substituent_pattern = Chem.MolFromSmarts('[cH]-[OH]')
    methoxy_pattern = Chem.MolFromSmarts('[cH]-O-C')
    has_oxy_substituent = mol.HasSubstructMatch(oxy_substituent_pattern)
    has_methoxy = mol.HasSubstructMatch(methoxy_pattern)
    if not (has_oxy_substituent or has_methoxy):
        return False, "No oxygenated substituents found on aromatic rings"

    # Check for positively charged oxygen atom
    has_positive_oxygen = any(atom.GetSymbol() == 'O' and atom.GetFormalCharge() == 1 for atom in mol.GetAtoms())
    if not has_positive_oxygen:
        return False, "No positively charged oxygen atom found"

    # Passed all checks
    return True, "Molecule is an anthocyanidin cation with flavylium core and oxygenated substituents"