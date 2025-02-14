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

    # Check if it's a simple atom (e.g., noble gases)
    if num_atoms == 1:
        atomic_num = mol.GetAtomWithIdx(0).GetAtomicNum()
        if atomic_num in [2, 10, 18, 36, 54, 86]: # Noble gases: He, Ne, Ar, Kr, Xe, Rn
             return True, "Noble gas atom"
        if atomic_num == 1: #Hydrogen
            return True, "Hydrogen atom"
        return False, "Single atom but not a noble gas."


    # Check for diatomic molecules
    if num_atoms == 2:
        atomic_num1 = mol.GetAtomWithIdx(0).GetAtomicNum()
        atomic_num2 = mol.GetAtomWithIdx(1).GetAtomicNum()

        if (atomic_num1 == 1 and atomic_num2 == 1) or \
           (atomic_num1 == 8 and atomic_num2 == 8) or \
           (atomic_num1 == 9 and atomic_num2 == 9) or \
           (atomic_num1 == 17 and atomic_num2 == 17): # H2, O2, F2, Cl2,
                return True, "Diatomic molecule of common gas-phase elements"
        if (atomic_num1 == 1 and atomic_num2 == 17) or \
            (atomic_num1 == 6 and atomic_num2 == 8) or \
             (atomic_num1 == 17 and atomic_num2 == 1): # HCl, CO, ClH
            return True, "Diatomic molecule of common gas-phase elements"



    # Small size rule: if too large/heavy not a gas
    if num_heavy_atoms > 6:
         return False, "Too many heavy atoms"

    # Check for common gas-phase elements (H, C, N, O, F, Cl, noble gases)
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in [1, 6, 7, 8, 9, 17, 2, 10, 18, 36, 54, 86]:
            return False, "Contains elements not commonly in gases"

    # Check for heavy groups, polar bonds or cyclic structures
    has_polar_bond = mol.HasSubstructMatch(Chem.MolFromSmarts('[C](=[O])')) or \
                        mol.HasSubstructMatch(Chem.MolFromSmarts('[C]-[F]')) or \
                        mol.HasSubstructMatch(Chem.MolFromSmarts('[C]-[Cl]')) or \
                        mol.HasSubstructMatch(Chem.MolFromSmarts('[O]-[H]')) or \
                        mol.HasSubstructMatch(Chem.MolFromSmarts('[N]-[H]'))

    if has_polar_bond:
         return False, "Contains polar bonds"


    cyclic_pattern = Chem.MolFromSmarts('[*1]~[*~*]~[*~*1]')
    if mol.HasSubstructMatch(cyclic_pattern):
            return False, "Contains cyclic structure"


    # Very short chains and linear molecules
    if num_rotatable_bonds > 3:
        return False, "Too many rotatable bonds"
    
    if mol_wt > 100:
        return False, "Molecular weight too high"

    return True, "Small molecule with common gas-phase elements and simple structure"