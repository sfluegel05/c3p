"""
Classifies: CHEBI:138675 gas molecular entity
"""
# gas_molecular_entity.py
"""
Classifies: CHEBI gas molecular entity
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_gas_molecular_entity(smiles: str):
    """
    Determines if a molecule is a main group molecular entity that is gaseous at STP.
    
    A gas molecular entity is any main group molecule that is gaseous at standard temperature and pressure (0°C and 100 kPa).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is a gas molecular entity, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define main group elements (atomic numbers)
    main_group_atomic_numbers = [
        1,  2,  3,  4,  5,  6,  7,  8,  9, 10,
       11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
       31, 32, 33, 34, 35, 36, 37, 38, 49, 50,
       51, 52, 53, 54, 55, 56, 81, 82, 83, 84,
       85, 86, 87, 88
    ]
    
    # Check if all atoms are main group elements
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in main_group_atomic_numbers:
            return False, f"Contains non-main group element: {atom.GetSymbol()}"
    
    # Calculate molecular weight
    mol_weight = Descriptors.MolWt(mol)
    
    # Estimate if the molecule is gaseous at STP
    # Using a heuristic that molecules with molecular weight less than 80 g/mol are likely gases
    if mol_weight <= 80:
        return True, f"Molecular weight {mol_weight:.2f} g/mol suggests the molecule is gaseous at STP"
    else:
        return False, f"Molecular weight {mol_weight:.2f} g/mol suggests the molecule is not gaseous at STP"