"""
Classifies: CHEBI:13601 3-oxo-5alpha-steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_oxo_5alpha_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5alpha-steroid based on its SMILES string.
    A 3-oxo-5alpha-steroid has a steroid skeleton, a carbonyl at position 3 and alpha
    configuration at position 5

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-5alpha-steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Define steroid core SMARTS pattern.
    # The numbers in the comments indicate the atom index based on the SMARTS matching
    # The order of the atoms is crucial for getting the right indices
    steroid_core_smarts = "[C]1[C]([H])([H])[C]([H])([H])[C]2([H])[C]([H])([H])[C]3([H])[C]([H])([H])[C]4([H])[C]([H])([H])[C]([H])([H])[C]([H])([H])[C]3([H])[C]5([H])[C]([H])([H])[C]([H])([H])[C]6([H])[C]([H])([H])[C]([H])([H])[C]([H])([H])[C]5([H])[C]1([H])([H])[C]([H])([H])[C]2([H])([H])[C]4([H])([H])6"
    
    steroid_core = Chem.MolFromSmarts(steroid_core_smarts)
    if steroid_core is None:
       return False, "Invalid core pattern"

    core_matches = mol.GetSubstructMatch(steroid_core)
    if not core_matches:
        return False, "Not a steroid: core structure not found"

    # Assign atom indices based on matching order.
    # We use the index values for carbon numbers 3 and 5 in the core
    # as they are assigned by the SMARTS pattern.
    atom_index_3 = core_matches[6] # index of carbon 3
    atom_index_5 = core_matches[13] # index of carbon 5


    # 2. Check for carbonyl at position 3
    atom3 = mol.GetAtomWithIdx(atom_index_3)
    if atom3.GetDegree() != 3:
        return False, "Carbon at position 3 not part of carbonyl"
    
    is_carbonyl_carbon_3 = False
    for neighbor in atom3.GetNeighbors():
      if neighbor.GetSymbol() == 'O' and mol.GetBondBetweenAtoms(atom3.GetIdx(), neighbor.GetIdx()).GetBondType() == Chem.rdchem.BondType.DOUBLE:
          is_carbonyl_carbon_3 = True
          break

    if not is_carbonyl_carbon_3:
        return False, "No carbonyl group at position 3 found."

    # 3. Check for alpha configuration at position 5
    atom5 = mol.GetAtomWithIdx(atom_index_5)

    chiral_tag = atom5.GetChiralTag()

    if chiral_tag != Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW and chiral_tag != Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW:
         return False, "No 5 alpha configuration found. Not a chiral carbon."

    # alpha configuration corresponds to CW or CCW, depending on the other bonds but using the RDKit standard
    # we can assume CW corresponds to alpha

    #  We need to check if the 5 atom has an explicit hydrogen and if it has chiral tag of CHI_TETRAHEDRAL_CW
    #  or not an explicit hydrogen and has chiral tag CHI_TETRAHEDRAL_CCW
    #  If these conditions are met we consider the 5 position to have the alpha configuration.
    has_explicit_hydrogen = False
    for neighbor in atom5.GetNeighbors():
        if neighbor.GetAtomicNum() == 1:
            has_explicit_hydrogen = True

    if has_explicit_hydrogen and chiral_tag == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW :
        pass
    elif not has_explicit_hydrogen and chiral_tag == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW:
        pass
    else:
       return False, "No 5 alpha configuration found. Wrong chirality or missing Hydrogen"

    return True, "Molecule is a 3-oxo-5alpha-steroid"