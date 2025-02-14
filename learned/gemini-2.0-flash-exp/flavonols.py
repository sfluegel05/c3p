"""
Classifies: CHEBI:28802 flavonols
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_flavonols(smiles: str):
    """
    Determines if a molecule is a flavonol based on its SMILES string.
    A flavonol is a flavone with a hydroxyl group at the 3-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the flavone core pattern using SMARTS, allowing for substitutions on the rings.
    #  'c1ccccc1-c2oc(-c3ccccc3)c(=O)c2' but with substitutions allowed
    flavone_core = Chem.MolFromSmarts('c1cc[c,C]c[c,C]c1-c2oc(-[c,C]3cc[c,C]c[c,C]c3)c(=O)c2')
    if not mol.HasSubstructMatch(flavone_core):
        return False, "Not a flavone core structure"

    # Check for the hydroxyl group at position 3 of the pyran ring (C ring).
    # This is the carbon directly attached to the carbonyl group in the flavone core
    # We first find the C=O group then check for the -C(OH)
    
    carbonyl_pattern = Chem.MolFromSmarts('C=O')
    carbonyl_matches = mol.GetSubstructMatches(carbonyl_pattern)
    found_oh_at_3 = False
    for match in carbonyl_matches:
        carbonyl_atom_idx = match[0]
        #get neighbors and find a carbon
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_atom_idx)
        for neighbor in carbonyl_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6: #Found a Carbon
                #check the next neighbor
                for neighbor2 in neighbor.GetNeighbors():
                    if neighbor2.GetAtomicNum() == 8 and neighbor2.GetFormalCharge() == 0 and neighbor2.GetTotalNumHs() ==1: #Found an oxygen and not a carbonyl
                        found_oh_at_3 = True
                        break
                if found_oh_at_3:
                    break #if found then we're done
        if found_oh_at_3:
            break
    
    if not found_oh_at_3:
      return False, "No hydroxyl group at position 3 of the flavone core"


    # If both the flavone core and 3-OH are found, it's a flavonol
    return True, "Flavonol: Flavone core with hydroxyl at position 3"