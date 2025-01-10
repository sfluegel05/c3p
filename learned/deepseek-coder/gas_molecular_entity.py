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
    main_group_elements = {1, 2, 5, 6, 7, 8, 9, 10, 13, 14, 15, 16, 17, 18, 36, 54, 86}  # Includes Kr, Xe, Rn
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in main_group_elements:
            return False, "Contains non-main group elements"

    # Calculate molecular weight
    mol_wt = Descriptors.MolWt(mol)

    # Check for noble gases (including isotopes)
    noble_gas_patterns = [
        "[He]", "[Ne]", "[Ar]", "[Kr]", "[Xe]", "[Rn]",
        "[*He]", "[*Ne]", "[*Ar]", "[*Kr]", "[*Xe]", "[*Rn]"
    ]
    if any(pattern in smiles for pattern in noble_gas_patterns):
        return True, "Noble gas"

    # Check for diatomic molecules (including isotopes)
    diatomic_patterns = [
        "[H][H]", "[1H][1H]", "[3H][3H]", "N#N", "[O][O]", "FF", "ClCl", "[C-]#[O+]",
        "I[H]", "Cl[H]", "Br[H]", "F[H]", "O=O", "N=N"
    ]
    if any(pattern in smiles for pattern in diatomic_patterns):
        return True, "Diatomic molecule"

    # Check for small molecules and common gases
    common_gases = {
        "O=C=O", "[H]N([H])[H]", "[O-][O+]=O", "C1CO1",  # CO2, NH3, O3, oxirane
        "C#C", "C=C", "C-C",  # Acetylene, ethene, ethane
        "C(=O)=O", "C#N", "C=N"  # CO2, HCN, etc.
    }
    if smiles in common_gases:
        return True, "Common gas molecule"

    # Check for small hydrocarbons (up to 4 carbons)
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[C,c]")):
        c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        if c_count <= 4 and mol_wt < 150:
            return True, "Small hydrocarbon"

    # Check for other small molecules (MW < 150 g/mol)
    if mol_wt < 150:
        return True, "Small molecule likely to be gas at STP"

    # If none of the above conditions are met, return None
    return None, "Cannot determine if the molecule is a gas at STP based on SMILES alone"