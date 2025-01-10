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
            return True, "Carbon atom"
            
    # Check molecular weight - most gases are light
    if mol_wt > 150:  # Most common gases are below this weight
        # Special case for perfluorinated compounds which can be gases at higher weights
        num_F = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[F]')))
        if num_F < 6:  # Allow higher weight only for heavily fluorinated compounds
            return False, f"Molecular weight too high ({mol_wt:.1f} Da)"
    
    # Most gases have relatively few heavy atoms
    if num_heavy_atoms > 8 and 'F' not in smiles:
        return False, f"Too many heavy atoms ({num_heavy_atoms})"
        
    # Common gas patterns
    gas_patterns = {
        'diatomic': '[!C!H]~[!C!H]',  # For N2, O2, F2, Cl2, etc.
        'carbon_oxide': '[C]=[O,N]',  # For CO, CO2
        'small_alkane': '[CX4]~[CX4]',  # For CH4, C2H6, etc.
        'small_alkene': '[CX3]=[CX3]',  # For C2H4, etc.
        'small_alkyne': '[CX2]#[CX2]',  # For C2H2, etc.
        'hydrogen_halide': '[F,Cl,Br,I][H]',  # For HF, HCl, HBr, HI
        'ammonia': '[NX3]([H])([H])[H]',  # NH3
        'ozone': '[O-][O+]=O'  # O3
    }
    
    for name, pattern in gas_patterns.items():
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            return True, f"Matches {name} pattern"
            
    # Check for small hydrocarbons (up to C4)
    if all(atom.GetSymbol() in ('C', 'H') for atom in mol.GetAtoms()):
        if num_heavy_atoms <= 4:
            return True, f"Small hydrocarbon with {num_heavy_atoms} carbons"
            
    # Check for simple molecules with few non-H atoms
    if num_heavy_atoms <= 3:
        return True, f"Simple molecule with {num_heavy_atoms} heavy atoms"
        
    # If we haven't returned True by now, it's probably not a gas
    return False, "Does not match known gas patterns and is too complex"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:33262',
        'name': 'gas molecular entity',
        'definition': 'Any main group molecular entity that is gaseous at standard temperature and pressure (STP; 0degreeC and 100 kPa).',
        'parents': []
    }
}