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
    
    # Check for single atoms first (noble gases, etc.)
    if num_heavy_atoms == 1:
        atom = mol.GetAtomWithIdx(0)
        symbol = atom.GetSymbol()
        
        # Noble gases (any isotope)
        if symbol in {'He', 'Ne', 'Ar', 'Kr', 'Xe', 'Rn'}:
            if atom.GetFormalCharge() != 0:
                return False, f"Charged {symbol} species"
            return True, f"Noble gas ({symbol})"
            
        # Single carbon atom (neutral only)
        if symbol == 'C' and atom.GetFormalCharge() == 0:
            return True, "Atomic carbon"
            
        return False, "Single atom but not a known gaseous element"

    # Check for ionic species (exclude most charged molecules)
    total_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if abs(total_charge) > 1:
        return False, "Multiple charges present"

    # Known diatomic gases and their isotopes
    diatomic_patterns = {
        '[H][H]': 'Hydrogen',
        'FF': 'Fluorine',
        'ClCl': 'Chlorine',
        '[O][O]': 'Oxygen',
        '[N]#[N]': 'Nitrogen'
    }
    
    for pattern, name in diatomic_patterns.items():
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            if num_heavy_atoms == 2 or (name == 'Hydrogen' and num_atoms == 2):
                return True, f"Diatomic gas: {name}"

    # Common small molecule gases
    if mol_wt <= 150:  # Increased threshold but with stricter checks
        # Simple inorganic gases
        inorganic_patterns = {
            'O=C=O': 'Carbon dioxide',
            '[C-]#[O+]': 'Carbon monoxide',
            '[O-][O+]=O': 'Ozone',
            '[H]N([H])[H]': 'Ammonia',
            'C1CO1': 'Oxirane',
            'F[H]': 'Hydrogen fluoride',
            'Cl[H]': 'Hydrogen chloride',
            'Br[H]': 'Hydrogen bromide',
            'I[H]': 'Hydrogen iodide'
        }
        
        for pattern, name in inorganic_patterns.items():
            if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)) and num_heavy_atoms <= 4:
                return True, f"Common inorganic gas: {name}"

        # Small hydrocarbons (C1-C4)
        if rdMolDescriptors.CalcNumLipinskiHBA(mol) <= 2 and \
           rdMolDescriptors.CalcNumLipinskiHBD(mol) <= 1:
            
            # Count carbon atoms
            carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
            
            if carbon_count <= 4:
                # Alkanes
                if '[C]' in Chem.MolToSmarts(mol) and not mol.HasSubstructMatch(Chem.MolFromSmarts('C=C')):
                    return True, f"Small alkane (C{carbon_count})"
                
                # Alkenes
                if mol.HasSubstructMatch(Chem.MolFromSmarts('C=C')) and not mol.HasSubstructMatch(Chem.MolFromSmarts('C(=O)')):
                    return True, f"Small alkene (C{carbon_count})"
                
                # Alkynes
                if mol.HasSubstructMatch(Chem.MolFromSmarts('C#C')):
                    return True, f"Small alkyne (C{carbon_count})"

    # Exclude if contains typical non-gas features
    non_gas_patterns = [
        '[OH]', # Hydroxyl groups (except in specific known gases)
        '[NH2]', # Primary amines (except ammonia)
        '[#6]~[#6]~[#6]~[#6]~[#6]', # Long carbon chains
        '[#6]1[#6][#6][#6][#6]1', # 5+ membered rings
        '[#6](=O)[OH]', # Carboxylic acids
        '[#6](=O)[O-]', # Carboxylates
        '[NH3+]', # Ammonium groups
        '[#15,#16]', # P, S (except in specific known gases)
    ]
    
    for pattern in non_gas_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            return False, "Contains functional groups typical of non-gases"

    return False, "Does not match known gas patterns or has characteristics of non-gases"