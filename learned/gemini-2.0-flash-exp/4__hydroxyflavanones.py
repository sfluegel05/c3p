"""
Classifies: CHEBI:140331 4'-hydroxyflavanones
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_4__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 4'-hydroxyflavanone based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
         bool: True if molecule is a 4'-hydroxyflavanone, False otherwise.
         str:  Reason for the classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for the flavanone core
    flavanone_core_pattern = Chem.MolFromSmarts("c1ccccc2C(=O)CC(Oc2c1)c1ccccc1")
    if not mol.HasSubstructMatch(flavanone_core_pattern):
        return False, "Molecule does not contain the flavanone core."
    
    # SMARTS pattern for the 4'-hydroxy substitution - O must be attached to C, and optionally H.
    hydroxy_4_pattern = Chem.MolFromSmarts("c1ccc(O)cc1")

    # Get the atoms matched by the flavanone core
    core_matches = mol.GetSubstructMatches(flavanone_core_pattern)
    
    # Iterate over the matched cores.
    for match in core_matches:
        # numbering of the core is: c1-c2-c3-c4-c5-c6-C7(=O)-C8-C9(Oc2-c1)-c10-c11-c12-c13-c14-c1
        # The atom at 4' position is the one at position 12
        
        # Get the 12th atom in the match (the one in position 4').
        atom_4prime_index = match[11]
        atom_4prime = mol.GetAtomWithIdx(atom_4prime_index)

        # Check if the atom has a direct Oxygen substituent
        for neighbor in atom_4prime.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:
                # Check if the oxygen has only one neighbor (i.e. is a hydroxyl)
                if len(neighbor.GetNeighbors()) == 1 :
                     
                    # If the oxygen is a hydroxyl, this molecule is a valid 4'-hydroxyflavanone
                    return True, "Molecule contains a flavanone core with a hydroxy substituent at the 4' position of the phenyl ring"

    # If no match is found, return false.
    return False, "Molecule does not have a hydroxyl substituent at the 4' position of the phenyl ring"