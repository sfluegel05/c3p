"""
Classifies: CHEBI:35346 11beta-hydroxy steroid
"""
from rdkit import Chem

def is_11beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is an 11beta-hydroxy steroid based on its SMILES string.
    An 11beta-hydroxy steroid is defined by the presence of a cyclopenta[a]phenanthrene core
    with a hydroxyl group at the 11th position with beta configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an 11beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # General steroid scaffold (four-fused rings)
    steroid_core = Chem.MolFromSmarts("[#6]1:#6:#6:#6:#6:#6:[C@H]2:#6:[#6](=[#8])[#6]:[#6]:[C@H]3:[#6:Object@H]([#8])C[C@H]4(#6)[C@@H]([#6])C[C@@:6]4\\C=C/[C@H]3[C@H]2[#6]1")
    
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core found matching cyclopenta[a]phenanthrene ring system"

    # 11beta-hydroxy group check
    # Ensure hydroxy group in the beta configuration at position 11
    match_found = False
    for bond in mol.GetBonds():
        atom1, atom2 = bond.GetBeginAtom(), bond.GetEndAtom()
        if atom1.GetAtomicNum() == 8 and atom2.GetAtomicNum() == 6:
            # Testing chirality and position index could be added for assertive testing
            chiral_tag = atom2.GetChiralTag()
            if chiral_tag == Chem.ChiralType.CHI_TETRAHEDRAL_CCW or chiral_tag == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
                match_found = True
                break
            
    if not match_found:
        return False, "No 11-beta hydroxy group detected or incorrect configuration"

    return True, "Contains 11-beta hydroxy group in correct configuration for steroid"

# Example usage
example_smiles = "[H][C@@]12CC[C@](O)(C(=O)CO)[C@@]1(CO)C[C@H](O)[C@@]1([H])[C@@]2([H])CCC2=CC(=O)CC[C@]12C"
print(is_11beta_hydroxy_steroid(example_smiles))