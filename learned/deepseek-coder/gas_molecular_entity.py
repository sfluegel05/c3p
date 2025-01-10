"""
Classifies: CHEBI:138675 gas molecular entity
"""
"""
Classifies: CHEBI:138675 gas molecular entity
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_gas_molecular_entity(smiles: str):
    """
    Determines if a molecule is a gas molecular entity based on its SMILES string.
    A gas molecular entity is any main group molecular entity that is gaseous at standard temperature and pressure (STP; 0°C and 100 kPa).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely a gas molecular entity, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule is a main group molecular entity
    # Main group elements are those in groups 1, 2, and 13-18 of the periodic table
    main_group_elements = {1, 2, 5, 6, 7, 8, 9, 10, 13, 14, 15, 16, 17, 18}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in main_group_elements:
            return False, "Contains non-main group elements"

    # Calculate molecular weight
    mol_wt = Descriptors.MolWt(mol)

    # Heuristic: Molecules with molecular weight < 100 g/mol are more likely to be gases at STP
    if mol_wt > 100:
        return False, f"Molecular weight {mol_wt:.2f} g/mol is too high for a gas at STP"

    # Check for common gas molecules (e.g., noble gases, diatomic molecules, small hydrocarbons)
    # Noble gases: He, Ne, Ar, Kr, Xe, Rn
    noble_gases = {"[He]", "[Ne]", "[Ar]", "[Kr]", "[Xe]", "[Rn]"}
    if smiles in noble_gases:
        return True, "Noble gas"

    # Diatomic molecules: H2, N2, O2, F2, Cl2
    diatomic_molecules = {"[H][H]", "N#N", "[O][O]", "FF", "ClCl"}
    if smiles in diatomic_molecules:
        return True, "Diatomic molecule"

    # Small hydrocarbons: methane, ethane, ethene, ethyne, propane, propene, etc.
    hydrocarbon_pattern = Chem.MolFromSmarts("[C,c]")
    if mol.HasSubstructMatch(hydrocarbon_pattern):
        # Check if the number of carbons is small (<= 4)
        c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        if c_count <= 4:
            return True, "Small hydrocarbon"

    # Other common gases: CO2, CO, NH3, etc.
    common_gases = {"O=C=O", "[C-]#[O+]", "[H]N([H])[H]"}
    if smiles in common_gases:
        return True, "Common gas molecule"

    # If none of the above conditions are met, return None
    return None, "Cannot determine if the molecule is a gas at STP based on SMILES alone"