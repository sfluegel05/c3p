"""
Classifies: CHEBI:16366 anthocyanidin cation
"""
"""
Classifies: CHEBI:36334 anthocyanidin cation
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_anthocyanidin_cation(smiles: str):
    """
    Determines if a molecule is an anthocyanidin cation based on its SMILES string.
    An anthocyanidin cation is an oxygenated derivative of flavylium (2-phenylchromenylium) with a positive charge,
    and is an aglycon (non-sugar part) of an anthocyanin cation.

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
    
    # Look for flavylium core (chromen-4-ylium)
    flavylium_pattern = Chem.MolFromSmarts("[o+]1c2ccccc2cc1")
    if not mol.HasSubstructMatch(flavylium_pattern):
        return False, "No flavylium core found"
    
    # Look for phenyl ring
    phenyl_pattern = Chem.MolFromSmarts("c1ccccc1")
    if not mol.HasSubstructMatch(phenyl_pattern):
        return False, "No phenyl ring found"
    
    # Check for positive charge
    formal_charge = rdMolDescriptors.CalcFormalCharge(mol)
    if formal_charge != 1:
        return False, "Incorrect formal charge (expected +1)"
    
    # Check for absence of sugar moieties
    sugar_pattern = Chem.MolFromSmarts("OC")
    if mol.HasSubstructMatch(sugar_pattern):
        return False, "Contains sugar moieties (not an aglycon)"
    
    # Look for oxygens attached to the flavylium core
    flavylium_atoms = mol.GetSubstructMatches(flavylium_pattern)[0]
    oxygens_on_core = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and any(bond.GetBeginAtomIdx() in flavylium_atoms or bond.GetEndAtomIdx() in flavylium_atoms for bond in atom.GetBonds())]
    if not oxygens_on_core:
        return False, "No oxygens attached to flavylium core"
    
    return True, "Contains oxygenated flavylium core with phenyl ring and positive charge, without sugar moieties"