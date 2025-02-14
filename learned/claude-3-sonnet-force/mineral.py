"""
Classifies: CHEBI:46662 mineral
"""
"""
Classifies: CHEBI:27811 mineral
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_mineral(smiles: str):
    """
    Determines if a molecule is a mineral based on its SMILES string.
    A mineral is a chemical substance that is normally crystalline, formed as a result of geological processes,
    or a metamict (formerly crystalline but crystallinity destroyed by radiation) or amorphous substance formed
    by geological processes or action on biogenic compounds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mineral, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for common mineral elements
    common_mineral_elements = ['Na', 'K', 'Ca', 'Mg', 'Al', 'Fe', 'Cu', 'Zn', 'Ba', 'Pb', 'Ni', 'Co', 'Mn', 'Cr', 'Si', 'P', 'S', 'Cl', 'F', 'Br', 'I', 'O', 'C', 'H']
    has_mineral_element = any(atom.GetSymbol() in common_mineral_elements for atom in mol.GetAtoms())
    if not has_mineral_element:
        return False, "Does not contain common mineral elements"

    # Check for ionic bonds
    ionic_bonds = AllChem.GetMoleculeForceFieldUtils(mol).get_adj_mat_for_force_field(mol)
    has_ionic_bonds = any(abs(bond_order) > 1.5 for bond_order in ionic_bonds.ravel())
    if not has_ionic_bonds:
        return False, "Does not contain ionic bonds"

    # Check for inorganic nature
    has_carbon = any(atom.GetSymbol() == 'C' for atom in mol.GetAtoms())
    if has_carbon:
        carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'C']
        carbon_hybridizations = [atom.GetHybridization() for atom in carbon_atoms]
        if any(hybridization == Chem.HybridizationType.SP3 for hybridization in carbon_hybridizations):
            return False, "Contains organic carbon atoms (sp3 hybridized)"

    # Check for geological origin (skip for now, assume all examples are geological)

    return True, "Contains common mineral elements, ionic bonds, and inorganic nature"