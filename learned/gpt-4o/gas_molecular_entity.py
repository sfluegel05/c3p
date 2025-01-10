"""
Classifies: CHEBI:138675 gas molecular entity
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_gas_molecular_entity(smiles: str):
    """
    Determines if a molecule is a gas molecular entity based on its SMILES string.
    A gas molecular entity is any main group molecular entity that is gaseous at
    standard temperature and pressure (STP; 0°C and 100 kPa).

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
    
    # Check for single atom gases (noble gases, specific light atoms)
    single_atom_gases = {'He', 'Ne', 'Ar', 'Kr', 'Xe', 'Rn', 'H'}
    if mol.GetNumAtoms() == 1:
        atom_symbol = mol.GetAtomWithIdx(0).GetSymbol()
        if atom_symbol in single_atom_gases:
            return True, f"Single atom gas: {atom_symbol}"
        else:
            return False, f"Single atom not typically a gas at STP: {atom_symbol}"
    
    # Check for common known small gas molecules
    small_gas_molecules = {
        'O=C=O', '[H]N([H])[H]', 'CC', 'CCC', 'C=C', 'CCCC', '[H]C([H])([H])[H]', '[O][O]', '[H][H]',
        'Cl[H]', 'CC(C)=C', 'CCC=C', '[O-][O+]=O', '[H][H]', 'CH4', 'N2', 'Cl2', 'HCl', 'I[H]',
        '[C-]#[O+]', 'O=O', 'F2', 'H2O', '[He]', '[Ne]', '[Ar]', 'C#C', 'C#N'
    }
    
    if smiles in small_gas_molecules:
        return True, "Recognized as a known small gas molecule at STP"

    # Calculate exact molecular weight and check molecule complexity
    molecular_weight = rdMolDescriptors.CalcExactMolWt(mol)
    ring_info = mol.GetRingInfo()
    
    # Check simple linear or small cyclic compounds with weight <= 60 Da
    if molecular_weight <= 60 and not ring_info.NumRings() > 1:
        return True, "Molecular weight and simplicity suggest it is likely a gas at STP"
    
    return False, "Too large or complex or charged to be a gas at STP"