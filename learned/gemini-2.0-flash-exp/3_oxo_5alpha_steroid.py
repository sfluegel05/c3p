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
    steroid_core_smarts = "[C]12[C][C]3[C]([C]4[C]([C]1[C]5[C]2[C]34)CC5)[C]"
    steroid_core = Chem.MolFromSmarts(steroid_core_smarts)
    if steroid_core is None:
        return False, "Invalid core pattern"
    core_matches = mol.GetSubstructMatch(steroid_core)
    if not core_matches:
        return False, "Not a steroid: core structure not found"

    # 2. Identify C3 using SMARTS relative to the core.
    c3_smarts = "[C]12[C][C](=[O])[C]3[C]([C]4[C]([C]1[C]5[C]2[C]34)CC5)[C]"
    c3_match = mol.GetSubstructMatch(Chem.MolFromSmarts(c3_smarts))
    if not c3_match:
         return False, "No carbonyl group at position 3 found."
    
    atom_index_3 = c3_match[2] # index of carbon 3

    # 3. Identify C5 using SMARTS relative to the core.
    c5_smarts = "[C]12[C][C]3[C]([C]4[C]([C]1[C]([C]5)[C]2[C]34)CC5)[C]"
    c5_match = mol.GetSubstructMatch(Chem.MolFromSmarts(c5_smarts))

    if not c5_match:
        return False, "Could not find C5"
    
    atom_index_5 = c5_match[8] # index of carbon 5


    # 4. Identify C4 using SMARTS relative to the core
    c4_smarts = "[C]12[C][C]3[C]([C]([C]4)[C]([C]1[C]5[C]2[C]34)CC5)[C]"
    c4_match = mol.GetSubstructMatch(Chem.MolFromSmarts(c4_smarts))
    if not c4_match:
        return False, "Could not find C4"
    atom_index_4 = c4_match[3]

    # 5. Identify C6 using SMARTS relative to the core
    c6_smarts = "[C]12[C][C]3[C]([C]4[C]([C]1[C]5[C]2[C]34)C([C]6)C5)[C]"
    c6_match = mol.GetSubstructMatch(Chem.MolFromSmarts(c6_smarts))

    if not c6_match:
         return False, "Could not find C6"
    atom_index_6 = c6_match[8]

    # 6. Check for alpha configuration at position 5
    atom5 = mol.GetAtomWithIdx(atom_index_5)
    atom4 = mol.GetAtomWithIdx(atom_index_4)
    atom6 = mol.GetAtomWithIdx(atom_index_6)

    # Check that the C5 atom has the correct chirality for an alpha configuration.
    # Specifically C5 bonded to C4 and C6 should have a clockwise configuration
    if not (atom5.GetChiralTag() == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW or atom5.GetChiralTag() == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW) :
        return False, "No 5 alpha configuration found. C5 is not a chiral center"

    #If we made it to here, then it is chiral, and the bonds to c4 and c6 should be oriented in a CW direction.
    # The default direction in RDKit for alpha is CW

    if atom5.GetChiralTag() == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW:
        return False, "No 5 alpha configuration found. Wrong chirality at C5"

    return True, "Molecule is a 3-oxo-5alpha-steroid"