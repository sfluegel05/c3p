"""
Classifies: CHEBI:67197 endocannabinoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_endocannabinoid(smiles: str):
    """
    Determines if a molecule is likely an endocannabinoid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely an endocannabinoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for long carbon chain (fatty acid or similar)
    fatty_acid_pattern = Chem.MolFromSmarts("C[C,C=C]~[C,C=C]~[C,C=C]~[C,C=C]~[C,C=C]")
    if not mol.HasSubstructMatch(fatty_acid_pattern):
        return False, "No characteristic fatty acid chain found"

    fatty_chain_matches = mol.GetSubstructMatches(Chem.MolFromSmarts("C[C,C=C]~[C,C=C]~[C,C=C]~[C,C=C]~[C,C=C]"))
    if not fatty_chain_matches:
        return False, "No characteristic fatty acid chain found"
    fatty_chain_carbons = len(fatty_chain_matches[0])  # Get number of carbons of the first match.
    if fatty_chain_carbons < 8:
        return False, "Chain too short to be a fatty acid"


    # 2. Check for polar head group (N or O atom connected to a carbon)
    has_polar_head = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in [7, 8]: #check for O or N
            for neighbor in atom.GetNeighbors():
                 if neighbor.GetAtomicNum() == 6: # if connected to a carbon
                     has_polar_head = True
                     break
            if has_polar_head:
                break
    if not has_polar_head:
        return False, "No polar head group (N or O) found"

    # 3. Check for connection between polar head and fatty acid chain ( a carbon that is linked to N or O, and other C's).
    head_group_found = False
    for atom in mol.GetAtoms():
       if atom.GetAtomicNum() == 6: #checks for carbon
          neighboring_atoms = [neighbor.GetAtomicNum() for neighbor in atom.GetNeighbors()]
          if (8 in neighboring_atoms or 7 in neighboring_atoms) and (6 in neighboring_atoms and neighboring_atoms.count(6) >=1):
             head_group_found=True
             break
    if not head_group_found:
       return False, "No link between chain and polar head"

    return True, "Likely an endocannabinoid based on long chain, polar head group, and linker"