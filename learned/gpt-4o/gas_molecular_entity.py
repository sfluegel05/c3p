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

    # Check for single atom gases (including isotopes)
    single_atom_gases = {'He', 'Ne', 'Ar', 'Kr', 'Xe', 'Rn', 'H'}
    isotopic_single_atoms = {f'[{isotope}{atom}]' for isotope in range(1, 300) for atom in single_atom_gases}
    if mol.GetNumAtoms() == 1:
        atom_symbol = mol.GetAtomWithIdx(0).GetSymbol()
        if atom_symbol in single_atom_gases or smiles in isotopic_single_atoms:
            return True, f"Single atom gas: {atom_symbol}"
        else:
            return False, f"Single atom not typically a gas at STP: {atom_symbol}"

    # Known small molecules that are gases at STP
    known_small_gases = [
        'O=C=O', '[C]', '[C-]#[O+]', '[O-][O+]=O', 'N#N', 'O=O', 'F[F]', 'ClCl',
        '[H]O[H]', '[H][H]', 'I[H]', '[H]N([H])[H]', 'Cl[H]', '[H]C([H])([H])[H]',
        'CC', 'CCC', 'CC(C)C', 'CCCC', 'C=C', 'CC=C', '[H]\C(C)=C(/[H])C', 'C1CO1'
    ]

    if smiles in known_small_gases:
        return True, "Recognized as a known small gas molecule at STP"
    
    # Calculate physicochemical properties
    molecular_weight = rdMolDescriptors.CalcExactMolWt(mol)
    ring_info = mol.GetRingInfo()

    # Check for molecule charge
    total_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if total_charge != 0:
        return False, "Molecule has net charge, unlikely to be a gas at STP"

    # Analyze simplicity and exclude complex charged or large compounds
    if molecular_weight <= 60 and (not ring_info.NumRings() > 1):
        return True, "Molecular weight and simplicity suggest it is likely a gas at STP"

    return False, "Too large, complex, or charged to be a gas at STP"