"""
Classifies: CHEBI:138675 gas molecular entity
"""
from rdkit import Chem
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
    
    # Check if it's a single atom (e.g., noble gases, H, or isotopes)
    if num_atoms == 1:
        atomic_num = mol.GetAtomWithIdx(0).GetAtomicNum()
        isotope = mol.GetAtomWithIdx(0).GetMass()
        if atomic_num in [1, 2, 10, 18, 36, 54, 86] or (atomic_num == 1 and isotope in [1.0078, 2.0141, 3.016]): #Include isotopes of H
            return True, "Noble gas or hydrogen atom"
        return False, "Single atom but not a noble gas or hydrogen or isotope."
    
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
           (atomic_num1 == 6 and atomic_num2 == 8) or \
           (atomic_num1 == 1 and atomic_num2 == 8) or \
           (atomic_num1 == 1 and atomic_num2 == 9) or \
           (atomic_num1 == 1 and atomic_num2 == 16) or \
           (atomic_num1 == 7 and atomic_num2 == 7) or \
            (atomic_num1 == 6 and atomic_num2 == 6): #Adding more diatomic combinations
           return True, "Diatomic molecule of common gas-phase elements"
    
    # Handle simple hydrides
    if num_heavy_atoms == 1 and num_atoms > 1:
        atomic_num = mol.GetAtomWithIdx(0).GetAtomicNum()
        if atomic_num in [1, 6, 7, 8, 9, 16, 17]:  #Main group elements that form gaseous hydrides
           return True, "Simple hydride gas"

    # Check for common gas-phase elements (non-metals and noble gases)
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in [1, 6, 7, 8, 9, 16, 17, 35, 53, 2, 10, 18, 36, 54, 86]: #Non metals and noble gases
           return False, "Contains elements not commonly found in gases"

    # Allow for more heavy atoms/larger molecules, and higher molecular weights.
    if num_heavy_atoms > 20:
        return False, "Too many heavy atoms for a typical gas"
    if mol_wt > 500:
        return False, "Molecular weight too high for a typical gas"

    #Less stringent functional group check
    functional_group_patterns = [
        Chem.MolFromSmarts("C(=O)O"),    # Carboxylic acid
        Chem.MolFromSmarts("C(=O)N"),    # Amide
       # Chem.MolFromSmarts("CO"),       # Alcohol - Removing for now
        Chem.MolFromSmarts("N"),         # Amine
        Chem.MolFromSmarts("NO2"),      # Nitro
        Chem.MolFromSmarts("S"),         # Thiol and thioether
        
        ]
    for pattern in functional_group_patterns:
        if mol.HasSubstructMatch(pattern):
            return False, "Contains functional groups that favor condensed phases"

    return True, "Likely a gas based on size, elements, and structure"