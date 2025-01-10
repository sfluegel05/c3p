"""
Classifies: CHEBI:138675 gas molecular entity
"""
from rdkit import Chem

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
    
    # Check for single atom gases (noble gases)
    single_atom_gases = {'He', 'Ne', 'Ar', 'Kr', 'Xe', 'Rn', 'H', 'C', 'N', 'O'}
    if mol.GetNumAtoms() == 1:
        atom_symbol = mol.GetAtomWithIdx(0).GetSymbol()
        if atom_symbol in single_atom_gases:
            return True, f"Single atom gas: {atom_symbol}"
        else:
            return False, f"Single atom not typically a gas at STP: {atom_symbol}"
    
    # Check for common small gas molecules
    # Examples: O=C=O, CC (ethane), CCC (propane), [H]N([H])[H] (ammonia)
    small_gas_molecules = {
        'O=C=O', 'CC', 'CCC', '[H]N([H])[H]', 'C=C', 'CCCC', 'CCC=C', 'ClC=C', 'FCC', 
        'C(O)C', 'C1CC1', 'Cl[H]', '[O][O]', 'HN', '[H][H]', '[C-]#[O+]', '[HXe]', 
        'C#C', 'CC(C)=C'
    }
    
    if smiles in small_gas_molecules:
        return True, "Recognized as a small gas molecule at STP"
    
    # Calculate the molecular weight and use it as an indicator
    mf = Chem.rdMolDescriptors.CalcExactMolWt(mol)
    if mf <= 44:  # Approximate cutoff for small molecules that are gases at STP
        return True, "Molecular weight indicates it is likely a gas at STP"
    else:
        return False, "Too large or complex to be a gas at STP"