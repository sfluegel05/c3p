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
    # This SMARTS pattern defines the steroid skeleton
    steroid_core_smarts = "[C]1[C][C]2[C][C]3[C][C]4[C]([C]3[C]5[C]([C]6[C]([C]4[C]21)CC[C]6)CC5)CC4"

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

    # Check for the alpha configuration based on the surrounding atoms.
    # In steroids, the alpha configuration at C5 has the C6 atom "below" the plane.
    # The alpha bond to the C6 at position 5 should be CW
    
    # Get C4 and C6 atoms
    atom_index_4 = core_matches[12] #index of carbon 4
    atom_index_6 = core_matches[14] #index of carbon 6
    atom4 = mol.GetAtomWithIdx(atom_index_4)
    atom6 = mol.GetAtomWithIdx(atom_index_6)
    
    #Get the bonds
    bond_c5_c4 = mol.GetBondBetweenAtoms(atom5.GetIdx(), atom4.GetIdx())
    bond_c5_c6 = mol.GetBondBetweenAtoms(atom5.GetIdx(), atom6.GetIdx())

    # Check that the C5 atom has the correct chirality for an alpha configuration.
    # Specifically C5 bonded to C4 and C6 should have a clockwise configuration
    if not (atom5.GetChiralTag() == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW or atom5.GetChiralTag() == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW) :
        return False, "No 5 alpha configuration found. C5 is not a chiral center"
        
    #If we made it to here, then it is chiral, and the bonds to c4 and c6 should be oriented in a CW direction.
    # The default direction in RDKit for alpha is CW

    if atom5.GetChiralTag() == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW:
        return False, "No 5 alpha configuration found. Wrong chirality at C5"
   
    return True, "Molecule is a 3-oxo-5alpha-steroid"