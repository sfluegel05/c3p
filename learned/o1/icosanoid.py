"""
Classifies: CHEBI:23899 icosanoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_icosanoid(smiles: str):
    """
    Determines if a molecule is an icosanoid based on its SMILES string.
    An icosanoid is a signaling molecule arising from oxidation of C20 essential fatty acids
    such as EPA, AA, and DGLA.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an icosanoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count total number of carbons in the molecule
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 18 or c_count > 24:
        return False, f"Total carbon count ({c_count}) not in range for icosanoids (18-24)"
    
    # Look for characteristic icosanoid substructures
    # Prostaglandin core: cyclopentane ring with side chains
    prostaglandin_pattern = Chem.MolFromSmarts("C1CCC(C1)C(=O)")
    has_prostaglandin = mol.HasSubstructMatch(prostaglandin_pattern)
    
    # Leukotriene pattern: linear chain with conjugated double bonds and hydroxyl groups
    leukotriene_pattern = Chem.MolFromSmarts("C=C/C=C/C=C/C(=O)")
    has_leukotriene = mol.HasSubstructMatch(leukotriene_pattern)
    
    # Thromboxane core: oxane (six-membered ring with one oxygen) with side chains
    thromboxane_pattern = Chem.MolFromSmarts("C1COC(C1)C(=O)")
    has_thromboxane = mol.HasSubstructMatch(thromboxane_pattern)
    
    if not (has_prostaglandin or has_leukotriene or has_thromboxane):
        return False, "No characteristic icosanoid substructures found"
    
    # Look for oxidative functional groups at specific positions
    hydroxyl_groups = mol.GetSubstructMatches(Chem.MolFromSmarts("[C;!$(C=O)][OH]"))
    keto_groups = mol.GetSubstructMatches(Chem.MolFromSmarts("C(=O)[C;!$(C=O)]"))
    epoxy_groups = mol.GetSubstructMatches(Chem.MolFromSmarts("C1OC1"))
    oxidative_mods = len(hydroxyl_groups) + len(keto_groups) + len(epoxy_groups)
    
    if oxidative_mods < 1:
        return False, "No oxidative modifications characteristic of icosanoids found"
    
    # Check for multiple double bonds
    unsat_bonds = [bond for bond in mol.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE]
    if len(unsat_bonds) < 3:
        return False, f"Only {len(unsat_bonds)} double bonds found, less than expected for icosanoids"
    
    return True, "Molecule matches key features of icosanoids: C20 backbone, characteristic substructures, oxidative modifications"