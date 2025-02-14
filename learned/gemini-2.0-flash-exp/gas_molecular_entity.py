"""
Classifies: CHEBI:138675 gas molecular entity
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_gas_molecular_entity(smiles: str):
    """
    Determines if a molecule is likely to be a gas at STP based on its SMILES string.
    This is an approximation based on rules and heuristics.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is likely a gas, False otherwise
        str: Reason for the classification
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    num_atoms = mol.GetNumAtoms()
    num_heavy_atoms = mol.GetNumHeavyAtoms()
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    num_rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    num_rings = Chem.GetSSSR(mol)


    # Check if it's a single atom (e.g., noble gases, H, or isotopes)
    if num_atoms == 1:
        atomic_num = mol.GetAtomWithIdx(0).GetAtomicNum()
        if atomic_num in [1, 2, 10, 18, 36, 54, 86]:
            return True, "Noble gas or hydrogen atom"
        return False, "Single atom but not a noble gas or hydrogen."
    
    # Check for diatomic molecules
    if num_atoms == 2:
       atomic_num1 = mol.GetAtomWithIdx(0).GetAtomicNum()
       atomic_num2 = mol.GetAtomWithIdx(1).GetAtomicNum()
       if (atomic_num1 == 1 and atomic_num2 == 1) or \
           (atomic_num1 == 8 and atomic_num2 == 8) or \
           (atomic_num1 == 9 and atomic_num2 == 9) or \
           (atomic_num1 == 17 and atomic_num2 == 17) or \
           (atomic_num1 == 1 and atomic_num2 == 17) or \
           (atomic_num1 == 1 and atomic_num2 == 53) or \
           (atomic_num1 == 6 and atomic_num2 == 8):
           return True, "Diatomic molecule of common gas-phase elements"
    
    # Handle single heavy atom with explicit hydrogens
    if num_heavy_atoms == 1 and num_atoms > 1:
        atomic_num = mol.GetAtomWithIdx(0).GetAtomicNum()
        if atomic_num in [6,7,17]: #C,N,Cl can form simple hydrides
           return True, "Simple hydride gas"

    # Check for common gas-phase elements (non-metals)
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in [1, 6, 7, 8, 9, 16, 17, 35, 53, 2, 10, 18, 36, 54, 86]: #Non metals and noble gases
           return False, "Contains elements not commonly found in gases"


    # Check the size and complexity of the molecule
    if num_heavy_atoms > 18:
        return False, "Too many heavy atoms for a typical gas"
    if num_atoms > 30:
        return False, "Too many atoms for a typical gas"
    if mol_wt > 350:
         return False, "Molecular weight too high for a typical gas"

    # Check for functional groups that tend to form condensed phases
    functional_group_patterns = [
        Chem.MolFromSmarts("C(=O)O"),    # Carboxylic acid
        Chem.MolFromSmarts("C(=O)N"),    # Amide
        Chem.MolFromSmarts("CO"),       # Alcohol
        Chem.MolFromSmarts("N"),         # Amine
        Chem.MolFromSmarts("C(=O)C"),    # Ketone
        Chem.MolFromSmarts("C(=O)H"),    # Aldehyde
        Chem.MolFromSmarts("NO2"),      # Nitro
        Chem.MolFromSmarts("C(=O)O"),    # Ester
        Chem.MolFromSmarts("S"),         # Thiol and thioether
        Chem.MolFromSmarts("O"),    # Ether
        Chem.MolFromSmarts("C=[O]"), # Carbonyl
        
    ]
    for pattern in functional_group_patterns:
        if mol.HasSubstructMatch(pattern):
            return False, "Contains functional groups that favor condensed phases"

    # Check for very long chains
    if num_rotatable_bonds > 8:
        return False, "Too many rotatable bonds for a typical gas"
    
    # Check for rings
    if num_rings > 2:
      return False, "Too many rings for a typical gas"

    return True, "Likely a gas based on size, elements, and structure"