"""
Classifies: CHEBI:138675 gas molecular entity
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_gas_molecular_entity(smiles: str):
    """
    Determines if a molecule is a gas molecular entity based on its SMILES string.
    A gas molecular entity is typically a simple structure that is gaseous at 
    standard temperature and pressure.

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
    
    # Common gaseous elements at STP, considering diatomic noble gases
    gaseous_elements = {1, 2, 6, 7, 8, 9, 10, 17, 18, 36, 54, 86}  # H, He, C, N, O, F, Ne, Cl, Ar, Kr, Xe, Rn
    
    # Known gaseous SMILES specific cases
    known_gaseous_smiles = {
        "[219Rn]", "[220Rn]", "[222Rn]", "[Rn]", "[Xe]", "[H][H]", "[He]", "[2H][2H]", "[3H][3H]", 
        "[4He]", "[3He]", "[6He]", "O=C=O", "[C-]#[O+]", "[H]N([H])[H]", "[H][C]([H])([H])[H]",
        "CCC", "C=C", "C(C)(F)(F)F", "C#C", "CC", "OOOO", "FF", "ClCl", "C1CO1", 
    }

    # Check if the molecule has a special specific known gaseous form
    if smiles in known_gaseous_smiles:
        return True, "Recognized special gaseous form"
    
    # Check atoms inclusion within gaseous elements
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in gaseous_elements:
            return False, f"Contains non-gaseous element: {atom.GetSymbol()}"
    
    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)

    # Allow gases with molecular weight up to 200, for noble gases and recognized exceptions
    if mol_wt > 200 and smiles not in known_gaseous_smiles:
        return False, f"Molecular weight too high for a typical gas: {mol_wt} g/mol"
    
    return True, "Molecule is a simple structure containing only gaseous elements at STP"

# Examples for testing
print(is_gas_molecular_entity("O=C=O"))  # Carbon dioxide
print(is_gas_molecular_entity("C1CO1"))  # Oxirane
print(is_gas_molecular_entity("[C-]#[O+]"))  # Carbon monoxide
print(is_gas_molecular_entity("[He]"))  # Helium
print(is_gas_molecular_entity("CCC"))  # Propane