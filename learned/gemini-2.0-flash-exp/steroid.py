"""
Classifies: CHEBI:35341 steroid
"""
from rdkit import Chem
from rdkit.Chem import rdmolops
from rdkit.Chem import AllChem


def is_steroid(smiles: str):
    """
    Determines if a molecule is a steroid based on its SMILES string.
    A steroid is defined as a molecule with a cyclopenta[a]phenanthrene skeleton,
    partially or completely hydrogenated, with methyl groups at C-10 and C-13,
    and often an alkyl group at C-17.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a generalized steroid core pattern using SMARTS (four fused rings)
    steroid_core_smarts = "[C]12[C]3[C]4[C]1[C]5[C]2[C]3[C]45"
    steroid_core_pattern = Chem.MolFromSmarts(steroid_core_smarts)

    # Check if the molecule has the basic four fused ring system
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "No steroid tetracyclic ring system found"

    # Check for methyl groups attached to the ring system using SMARTS and get the atoms in the ring system
    match = mol.GetSubstructMatch(steroid_core_pattern)
    atoms = [mol.GetAtomWithIdx(x) for x in match]

    # Get the neighbors of each ring atom that are not in the ring system
    neighbors = []
    for i,atom in enumerate(atoms):
      for neighbor in atom.GetNeighbors():
         if neighbor.GetIdx() not in match:
            neighbors.append((i,neighbor))
    
    methyl_atoms = []
    for (i,atom) in neighbors:
      if atom.GetSymbol() == 'C' and len(atom.GetNeighbors())==1:
        methyl_atoms.append( (i,atom) )

    # Check for at least two methyl groups connected to the ring system
    if len(methyl_atoms) < 2:
        return False, f"Found {len(methyl_atoms)} methyl groups on the ring system, at least two are required."

    return True, "Contains the steroid tetracyclic ring system with at least two methyl groups"