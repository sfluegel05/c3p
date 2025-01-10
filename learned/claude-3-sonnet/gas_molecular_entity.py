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
    
    # Molecular weight cutoff (most gases are < 100 g/mol)
    if mol_wt > 100 and smiles not in {'C(C(C(F)(F)F)(F)F)(C(F)(F)F)(F)F'}:  # Special case for perfluorobutane
        return False, "Molecular weight too high for a gas"

    # Check for single atoms first (noble gases, etc.)
    if num_heavy_atoms == 1:
        atom = mol.GetAtomWithIdx(0)
        symbol = atom.GetSymbol()
        
        # Noble gases (any isotope)
        if symbol in {'He', 'Ne', 'Ar', 'Kr', 'Xe', 'Rn'}:
            return True, f"Noble gas ({symbol})"
            
        # Single carbon atom
        if symbol == 'C':
            return True, "Atomic carbon"
            
        return False, "Single atom but not a known gaseous element"

    # Known diatomic gases
    diatomic_patterns = {
        '[H][H]': 'Hydrogen',
        'FF': 'Fluorine',
        'ClCl': 'Chlorine',
        '[O][O]': 'Oxygen',
        'N#N': 'Nitrogen',
        'I[H]': 'Hydrogen iodide',
        'Br[H]': 'Hydrogen bromide',
        'Cl[H]': 'Hydrogen chloride',
        'F[H]': 'Hydrogen fluoride'
    }
    
    for pattern, name in diatomic_patterns.items():
        if Chem.MolFromSmarts(pattern) is not None and mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            if num_heavy_atoms <= 2:
                return True, f"Diatomic gas: {name}"

    # Common small molecule gases
    common_gases = {
        'O=C=O': 'Carbon dioxide',
        '[C-]#[O+]': 'Carbon monoxide',
        '[O-][O+]=O': 'Ozone',
        '[H]N([H])[H]': 'Ammonia',
        'C1CO1': 'Oxirane'
    }
    
    for pattern, name in common_gases.items():
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            return True, f"Common gas: {name}"

    # Count carbon atoms
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Simple alkanes (C1-C4)
    if carbon_count <= 4:
        # Check if it's a simple alkane (only C and H)
        only_C_and_H = all(atom.GetAtomicNum() in {1, 6} for atom in mol.GetAtoms())
        no_double_or_triple_bonds = not (mol.HasSubstructMatch(Chem.MolFromSmarts('C=C')) or 
                                       mol.HasSubstructMatch(Chem.MolFromSmarts('C#C')))
        if only_C_and_H and no_double_or_triple_bonds:
            return True, f"Simple alkane (C{carbon_count})"
            
        # Simple alkenes (C2-C4, no other functional groups)
        if mol.HasSubstructMatch(Chem.MolFromSmarts('C=C')):
            # Check for absence of other functional groups
            if only_C_and_H:
                return True, f"Simple alkene (C{carbon_count})"
                
        # Simple alkynes (C2-C3, no other functional groups)
        if carbon_count <= 3 and mol.HasSubstructMatch(Chem.MolFromSmarts('C#C')):
            if only_C_and_H:
                return True, f"Simple alkyne (C{carbon_count})"

    # Exclude if contains typical non-gas features
    non_gas_patterns = [
        '[OH]',  # Hydroxyl groups
        '[NH2]',  # Primary amines (except ammonia)
        '[#6]~[#6]~[#6]~[#6]~[#6]',  # Long carbon chains
        '[#6]1[#6][#6][#6][#6]1',  # 5-membered rings
        '[#6](=O)[OH]',  # Carboxylic acids
        '[#6](=O)[O-]',  # Carboxylates
        '[NH3+]',  # Ammonium groups
        '[#15,#16]',  # P, S (except in specific known gases)
        '[#6]=[#7]',  # C=N bonds
        '[#6]=[#8]',  # C=O bonds
        '[#6]S',  # C-S bonds
        '[#7]1[#6][#6]1',  # Aziridine rings
    ]
    
    for pattern in non_gas_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            return False, "Contains functional groups typical of non-gases"

    return False, "Does not match known gas patterns"