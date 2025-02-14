"""
Classifies: CHEBI:1624 3-oxo-5beta-steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_5beta_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5beta-steroid based on its SMILES string.
    
    A 3-oxo-5beta-steroid has:
        - a steroid core
        - a ketone at the 3 position
        - a beta configuration at position 5.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-5beta-steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for steroid core
    steroid_core_pattern = Chem.MolFromSmarts("C1CC2CCC3C4(C)C(C)CCC(C)C(C)=C3C(C)C2C1") # basic 4-ring structure
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "Molecule does not contain a steroid core"
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    if num_rings != 4:
         return False, "Molecule does not contain a four-ring core"

    # 2. Check for 3-oxo group (C=O at position 3)
    # We define an atom list to correspond to the carbon positions in the ring pattern:
    # 1 is at position 1 in the pattern above
    # 2 is at position 2 in the pattern above
    # ...
    # 10 is the bridgehead C connecting rings 5 and 6 in the pattern.
    # In the SMARTS definition above, this carbon is the one at position "C2"
    # We need the 3rd carbon in the sequence, which is "CC2".
    
    # Find atom indices matching the SMARTS and use the first one in the match
    core_match = mol.GetSubstructMatch(steroid_core_pattern)
    if not core_match:
        return False, "Steroid core not found"

    # Map the steroid core atom indices from the pattern to the current mol
    atom_index_3 = core_match[2]  # index 2 refers to the carbon at position 3 in the steroid_core_pattern

    
    # Check for ketone at position 3 (atom_index_3)
    atom_3 = mol.GetAtomWithIdx(atom_index_3)

    if atom_3.GetTotalValence() != 3:
        return False, "Carbon at position 3 is not sp2"

    has_oxo = False
    for neighbor in atom_3.GetNeighbors():
        if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalValence() == 2: #double bonded O
            has_oxo = True
            break
    
    if not has_oxo:
        return False, "No ketone at position 3"

    # 3. Check for 5-beta configuration.
    # Find the atom at position 5, corresponding to the atom index 5 in the original steroid_core_pattern
    atom_index_5 = core_match[4]

    atom_5 = mol.GetAtomWithIdx(atom_index_5)

    # Check if any of the neighbors is a hydrogen.
    for neighbor in atom_5.GetNeighbors():
        if neighbor.GetAtomicNum() == 1: # this is an H
          bond = mol.GetBondBetweenAtoms(atom_5.GetIdx(),neighbor.GetIdx())
          if bond.GetStereo() == Chem.BondStereo.STEREOUP: # This is beta
            return True, "Molecule is a 3-oxo-5beta-steroid"
          else:
            return False, "Hydrogen at position 5 is not beta"


    return False, "No hydrogen found at position 5"