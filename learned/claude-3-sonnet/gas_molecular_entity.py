"""
Classifies: CHEBI:138675 gas molecular entity
"""
"""
Classifies: CHEBI:33262 gas molecular entity
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_gas_molecular_entity(smiles: str):
    """
    Determines if a molecule is likely to be gaseous at STP (0°C and 100 kPa)
    based on its molecular structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely a gas at STP, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get basic molecular properties
    mol_wt = Descriptors.ExactMolWt(mol)
    num_atoms = mol.GetNumAtoms()
    num_heavy_atoms = rdMolDescriptors.CalcNumHeavyAtoms(mol)
    
    # Check for noble gases (He, Ne, Ar, Kr, Xe, Rn)
    noble_gases = {'He', 'Ne', 'Ar', 'Kr', 'Xe', 'Rn'}
    if num_atoms == 1:
        symbol = mol.GetAtomWithIdx(0).GetSymbol()
        if symbol in noble_gases:
            return True, f"Noble gas ({symbol})"
        if symbol == 'C':
            return True, "Atomic carbon"

    # Most gases have very low molecular weight
    if mol_wt > 100:  # Lowered from 150
        # Special case for perfluorinated compounds
        num_F = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[F]')))
        if num_F < 8:  # More strict fluorine requirement
            return False, f"Molecular weight too high ({mol_wt:.1f} Da)"

    # Exclude common functional groups that make molecules non-gaseous at STP
    non_gas_patterns = [
        '[OH]C(=O)[!$(*[OH,NH2])]',  # Carboxylic acids
        '[OH]C([!$(C=O)])',  # Alcohols (except formaldehyde)
        '[NH2][!$(C=O)][!$(C=O)]',  # Primary amines (except simple ones)
        '[#6]~[#6]~[#6]~[#6]~[#6]',  # Long carbon chains
        '[#6]1[#6][#6][#6][#6]1',  # Cyclopentane and larger rings
    ]
    
    for pattern in non_gas_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            return False, "Contains functional groups typical of non-gases"

    # Known gas patterns with specific constraints
    gas_patterns = {
        'diatomic': [
            '[H][H]',  # H2
            '[N]#[N]',  # N2
            '[O][O]',  # O2
            'FF',      # F2
            'ClCl',    # Cl2
            'II',      # I2
        ],
        'carbon_oxide': [
            'O=C=O',   # CO2
            '[C-]#[O+]',  # CO
        ],
        'small_alkane': [
            '[CH4]',     # Methane
            '[CH3][CH3]',  # Ethane
            '[CH3][CH2][CH3]',  # Propane
            '[CH3][CH2][CH2][CH3]',  # Butane
        ],
        'small_alkene': [
            'C=C',  # Ethene
            'CC=C',  # Propene
            'C/C=C/C',  # But-2-ene
        ],
        'small_alkyne': [
            'C#C',  # Ethyne
            'CC#C',  # Propyne
        ],
        'simple_molecules': [
            '[NH3]',  # Ammonia
            'C1CO1',  # Oxirane
            '[O-][O+]=O',  # Ozone
        ],
        'hydrogen_halides': [
            'F[H]',  # HF
            'Cl[H]',  # HCl
            'Br[H]',  # HBr
            'I[H]',   # HI
        ]
    }

    for category, patterns in gas_patterns.items():
        for pattern in patterns:
            pattern_mol = Chem.MolFromSmarts(pattern)
            if pattern_mol and mol.HasSubstructMatch(pattern_mol):
                # Additional check for size
                if num_heavy_atoms <= 5:
                    return True, f"Known gas: matches {category} pattern"

    # If we haven't returned True by now, it's probably not a gas
    return False, "Does not match known gas patterns or has characteristics of non-gases"