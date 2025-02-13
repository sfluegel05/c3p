"""
Classifies: CHEBI:138675 gas molecular entity
"""
"""
Classifies: CHEBI:51869 gas molecular entity
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_gas_molecular_entity(smiles: str):
    """
    Determines if a molecule is a gas molecular entity based on its SMILES string.
    Gas molecular entities are main group molecular entities that are gaseous at standard
    temperature and pressure (STP; 0°C and 100 kPa).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a gas molecular entity, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Exclude charged species
    if any(atom.GetFormalCharge() != 0 for atom in mol.GetAtoms()):
        return False, "Contains charged atoms, cannot be a gas molecular entity"
    
    # Exclude transition metals
    periodic_table = Chem.GetPeriodicTable()
    if any(periodic_table.GetNOuterElecs(atom.GetAtomicNum()) > 4 for atom in mol.GetAtoms()):
        return False, "Contains transition metals, cannot be a gas molecular entity"
    
    # Check molecular weight - gas molecular entities typically <200 Da
    mol_wt = Descriptors.MolWt(mol)
    if mol_wt > 200:
        return False, "Molecular weight too high for a gas molecular entity"
    
    # Check boiling point - gas molecular entities typically boil <100°C
    bp = Descriptors.MolWt(mol) * 0.57 + 198.2  # Approximate boiling point formula
    if bp > 373.15:  # 100°C in K
        return False, "Boiling point too high for a gas molecular entity"
    
    # Check for common gases
    smiles_gases = ['[He]', '[Ne]', '[Ar]', '[Kr]', '[Xe]', '[Rn]', '[H][H]', '[O][O]', 'O=C=O', 'C(C)(C)(C)', 'O=N=O', '[CH4]', '[NH3]', '[H][Cl]', '[H][F]', '[H][Br]', '[H][I]', '[C-]#[O+]', '[N-]#[N+]', '[F][F]', '[Cl][Cl]', '[Br][Br]', '[I][I]']
    if smiles in smiles_gases:
        return True, "Common gas molecular entity"
    
    # Check for alkanes, alkenes, alkynes, cycloalkanes, and ethers
    if all(atom.GetAtomicNum() in [1, 6, 8] for atom in mol.GetAtoms()):
        return True, "Alkane, alkene, alkyne, cycloalkane, or ether gas molecular entity"
    
    # If none of the above conditions are met, assume it's not a gas molecular entity
    return False, "Molecule does not appear to be a gas molecular entity"