"""
Classifies: CHEBI:138675 gas molecular entity
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_gas_molecular_entity(smiles: str):
    """
    Determines if a molecule is a gas molecular entity based on its SMILES string.
    A gas molecular entity is typically a simple structure that is gaseous at standard temperature and pressure.

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

    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    
    # Common gaseous elements
    gaseous_elements = {2, 6, 7, 8, 9, 10, 17, 18, 36, 54, 86} # H, C, N, O, F, Ne, Cl, Ar, Kr, Xe, Rn

    # Check atoms and molecular weight
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in gaseous_elements:
            return False, f"Contains non-gaseous element: {atom.GetSymbol()}"
    
    # Typical range for gaseous molecular weight
    if mol_wt > 100:
        return False, f"Molecular weight too high for a typical gas: {mol_wt} g/mol"
    
    return True, "Molecule is a simple structure containing only gaseous elements at STP"

# Examples for testing
print(is_gas_molecular_entity("O=C=O"))  # Carbon dioxide
print(is_gas_molecular_entity("C1CO1"))  # Oxirane
print(is_gas_molecular_entity("[C-]#[O+]"))  # Carbon monoxide
print(is_gas_molecular_entity("[He]"))  # Helium
print(is_gas_molecular_entity("CCC"))  # Propane