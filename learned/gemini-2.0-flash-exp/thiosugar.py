"""
Classifies: CHEBI:73754 thiosugar
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_thiosugar(smiles: str):
    """
    Determines if a molecule is a thiosugar based on its SMILES string.
    A thiosugar is a carbohydrate derivative in which one or more of the oxygens
    or hydroxy groups of the parent carbohydrate is replaced by sulfur or -SR,
    where R can be hydrogen or any group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a thiosugar, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of sulfur atoms
    sulfur_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16]
    if not sulfur_atoms:
        return False, "No sulfur atoms found"
    
    # Check for at least one oxygen atom
    oxygen_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8]
    if not oxygen_atoms:
        return False, "No oxygen atoms found (not a carbohydrate derivative)"

    # Detect ring systems
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()

    if not rings:
        # If no rings found, check for -SR type structures
        for s_atom in sulfur_atoms:
            for neighbor in s_atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:
                    return True, "Sulfur attached to a carbon, potentially -SR"
        return False, "No rings found, and sulfur not attached to C"

    # Check if any ring contains at least one C, one O and one S
    ring_with_S = False
    for ring in rings:
        has_carbon = False
        has_oxygen = False
        has_sulfur = False
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 6:
                has_carbon = True
            elif atom.GetAtomicNum() == 8:
                has_oxygen = True
            elif atom.GetAtomicNum() == 16:
                 has_sulfur = True
        if has_carbon and has_oxygen and has_sulfur:
            ring_with_S = True
            break

    if ring_with_S:
      return True, "Ring system contains C, O and S atoms"
    else:
        # Check if any sulfur atom is attached to a ring carbon
        for ring in rings:
           for s_atom in sulfur_atoms:
            for neighbor in s_atom.GetNeighbors():
               if neighbor.GetIdx() in ring:
                carbon_in_ring = False
                for ring_atom in ring:
                    atom = mol.GetAtomWithIdx(ring_atom)
                    if atom.GetAtomicNum() == 6 and neighbor.GetIdx() == atom.GetIdx():
                       carbon_in_ring = True
                       break
                if carbon_in_ring:
                  return True, "Sulfur attached to a carbon atom in a ring"
    return False, "No ring with C,O,S, or S attached to ring C"