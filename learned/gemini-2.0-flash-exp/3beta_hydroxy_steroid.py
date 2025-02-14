"""
Classifies: CHEBI:36836 3beta-hydroxy steroid
"""
"""
Classifies: 3beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy steroid based on its SMILES string.
    A 3beta-hydroxy steroid is a steroid with a hydroxyl group in the beta position on the 3rd carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more general SMARTS pattern for the steroid core with the 3-beta hydroxyl.
    # Here we look for 4 fused rings with a carbon atom at position 3 connected with a beta-bond to an oxygen.
    # We use an atom map so we can verify that we have the correct connectivity.
    steroid_pattern = Chem.MolFromSmarts("[C]12~[C]~[C]~[C]1~[C]~[C]3~[C]2~[C]~[C]4~[C]3~[C]([O;H1]@[H])@[C]4")

    if steroid_pattern is None:
      return False, "Invalid SMARTS pattern"

    # Check if the molecule matches the pattern
    matches = mol.GetSubstructMatches(steroid_pattern)

    if not matches:
        return False, "Molecule does not match the steroid core with beta-hydroxyl at position 3"
    
    for match in matches:
        # Get atom indices
        atom_indices = list(match)
        
        # Extract the atom mapped to the beta-hydroxyl substituent
        oxy_atom_index = atom_indices[-1]
        oxy_atom = mol.GetAtomWithIdx(oxy_atom_index)

        # Check that the oxygen is connected to the carbon atom specified by the substructure match
        # and check that is indeed has beta stereochemistry
        carbon_atom_index = atom_indices[-2]
        carbon_atom = mol.GetAtomWithIdx(carbon_atom_index)
        is_beta = False
        for bond in carbon_atom.GetBonds():
            if bond.GetOtherAtomIdx(carbon_atom.GetIdx()) == oxy_atom.GetIdx():
                if bond.GetBondDir() == Chem.rdchem.BondDir.BEGINWEDGE:
                  is_beta = True
                  break
        if not is_beta:
           continue # Check other matches

        # If we've reached here, the molecule matches the criteria
        return True, "Molecule matches the criteria for a 3beta-hydroxy steroid"

    return False, "Molecule does not match the steroid core with beta-hydroxyl at position 3"