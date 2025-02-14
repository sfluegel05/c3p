"""
Classifies: CHEBI:36141 quinone
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_quinone(smiles: str):
    """
    Determines if a molecule is a quinone based on its SMILES string.
    A quinone is a cyclic dione structure derived from an aromatic compound, with conjugated double bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quinone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get all rings
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()

    if not rings:
       return False, "No rings found"

    # Check if molecule has any aromatic rings
    has_aromatic = any(ring_info.IsRingAromatic(ring_idx) for ring_idx in range(len(rings)))
    if not has_aromatic:
       return False, "No aromatic rings found"

    # Define SMARTS pattern for quinone substructure
    # This pattern looks for a ring with at least two carbonyl groups (=O) separated by a conjugated double bond
    quinone_pattern1 = Chem.MolFromSmarts("[C:1](=[O:2])[C:3]=[C:4][C:5](=[O:6])")
    quinone_pattern2 = Chem.MolFromSmarts("[C:1](=[O:2])[C:3]=[C:4]~[C:5]=[C:6][C:7](=[O:8])")


    # Check each aromatic ring for the quinone pattern
    for ring in rings:
      if not ring_info.IsRingAromatic(ring):
        continue

      ring_atoms = [mol.GetAtomWithIdx(atom) for atom in ring]

      #create a submolecule only from the atoms of the current ring
      ring_mol = Chem.RWMol()
      for atom in ring_atoms:
          new_atom = Chem.Atom(atom.GetAtomicNum())
          ring_mol.AddAtom(new_atom)
      
      for i in range(len(ring_atoms)):
        for j in range(i + 1, len(ring_atoms)):
            if mol.GetBondBetweenAtoms(ring_atoms[i].GetIdx(), ring_atoms[j].GetIdx()):
               ring_mol.AddBond(i, j, mol.GetBondBetweenAtoms(ring_atoms[i].GetIdx(), ring_atoms[j].GetIdx()).GetBondType())

      ring_mol = ring_mol.GetMol()


      if ring_mol.HasSubstructMatch(quinone_pattern1) or ring_mol.HasSubstructMatch(quinone_pattern2):
          return True, "Aromatic ring with dione substructure detected"

    return False, "No quinone structure detected"