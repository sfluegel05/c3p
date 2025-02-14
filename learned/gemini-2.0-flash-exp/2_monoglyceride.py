"""
Classifies: CHEBI:17389 2-monoglyceride
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_2_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 2-monoglyceride based on its SMILES string.
    A 2-monoglyceride is a glycerol with one fatty acid attached to the central (2nd) carbon via an ester bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-monoglyceride, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Glycerol backbone pattern (C-C-C with oxygens)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4]([OX2])[CH2X4]")
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    if not glycerol_matches:
       return False, "No glycerol backbone found"
    
    # Check for *one* ester group
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
       return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Get the atoms of the glycerol backbone
    glycerol_atoms = list(glycerol_matches[0])
    
    # Check that the ester is on the middle C
    ester_connected = False
    ester_o_atom_index = ester_matches[0][0]
    ester_o_atom = mol.GetAtomWithIdx(ester_o_atom_index)
        
    for neighbor in ester_o_atom.GetNeighbors():
        if neighbor.GetIdx() == glycerol_atoms[1]:  # glycerol atoms: 0-1-2, so 1 is the middle carbon
            ester_connected = True
            break
    if not ester_connected:
        return False, "Ester group not connected to the central carbon of the glycerol"

    # Verify that the other 2 hydroxyls are free. Check that oxygens at C1 and C3 are part of hydroxyls
    free_OH_count = 0
    glycerol_carbons = [mol.GetAtomWithIdx(glycerol_atoms[0]), mol.GetAtomWithIdx(glycerol_atoms[2])] # carbons 1 and 3
    for carbon in glycerol_carbons:
        for neighbor in carbon.GetNeighbors():
            if neighbor.GetAtomicNum() == 8: # found an oxygen
                if len(neighbor.GetNeighbors()) == 1:  # This oxygen is bound to only one other atom (the carbon)
                    free_OH_count += 1
                else:
                     # if there are multiple neighbors, ensure that the oxygen not part of the ester group (connected to carbonyl)
                    is_ester_oxygen = False
                    for n2 in neighbor.GetNeighbors():
                       if n2.GetIdx() == ester_o_atom_index:
                         is_ester_oxygen = True
                         break
                    if not is_ester_oxygen:
                        if len(neighbor.GetNeighbors()) ==1:
                           free_OH_count += 1

    if free_OH_count != 2:
        return False, "Did not find two free hydroxyl groups at the 1 and 3 positions of the glycerol"


    return True, "Contains a glycerol backbone with one fatty acid attached at the 2nd carbon via ester bond"