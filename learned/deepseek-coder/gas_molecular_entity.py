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

    # Heuristic: Molecules with molecular weight < 200 g/mol are more likely to be gases at STP
    if mol_wt > 200:
        return False, f"Molecular weight {mol_wt:.2f} g/mol is too high for a gas at STP"

    # Check for noble gases (including isotopes)
    noble_gases = {"[He]", "[Ne]", "[Ar]", "[Kr]", "[Xe]", "[Rn]",
                  "[4He]", "[3He]", "[6He]", "[219Rn]", "[220Rn]", "[222Rn]"}
    if smiles in noble_gases:
        return True, "Noble gas"

    # Check for diatomic molecules (including isotopes)
    diatomic_molecules = {"[H][H]", "[1H][1H]", "[3H][3H]", "N#N", "[O][O]", "FF", "ClCl", "[C-]#[O+]"}
    if smiles in diatomic_molecules:
        return True, "Diatomic molecule"

    # Check for small hydrocarbons and common gases
    # Include molecules with up to 8 carbons and molecular weight < 200 g/mol
    hydrocarbon_pattern = Chem.MolFromSmarts("[C,c]")
    if mol.HasSubstructMatch(hydrocarbon_pattern):
        c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        if c_count <= 8 and mol_wt < 200:
            return True, "Small hydrocarbon"

    # Other common gases
    common_gases = {"O=C=O", "[H]N([H])[H]", "O", "N", "C", "I[H]", "Cl[H]"}
    if smiles in common_gases:
        return True, "Common gas molecule"

    # If none of the above conditions are met, return None
    return None, "Cannot determine if the molecule is a gas at STP based on SMILES alone"