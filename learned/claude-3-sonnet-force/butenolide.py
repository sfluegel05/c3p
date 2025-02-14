"""
Classifies: CHEBI:50523 butenolide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_butenolide(smiles: str):
    """
    Determines if a molecule is a butenolide based on its SMILES string.
    A butenolide is a gamma-lactone consisting of a 2-furanone skeleton and its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a butenolide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find furanone ring
    ring_info = mol.GetRingInfo()
    furanone_ring = None
    for ring in ring_info.AtomRings():
        if ring.Count() == 5:  # 5-membered ring
            atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
            has_carbonyl = any(atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == 0 and sum(bond.GetBondTypeAsDouble() for bond in atom.GetBonds()) == 2 for atom in atoms)
            has_oxygen = any(atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == 0 and sum(bond.GetBondTypeAsDouble() for bond in atom.GetBonds()) == 1 for atom in atoms)
            has_double_bond = any(bond.GetBondType() == Chem.BondType.DOUBLE for atom in atoms for bond in atom.GetBonds())
            if has_carbonyl and has_oxygen and has_double_bond:
                furanone_ring = ring
                break
    if furanone_ring is None:
        return False, "No furanone ring found"
    
    # Check gamma-lactone arrangement
    ring_atoms = [mol.GetAtomWithIdx(idx) for idx in furanone_ring]
    carbonyl_atom = next((atom for atom in ring_atoms if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == 0 and sum(bond.GetBondTypeAsDouble() for bond in atom.GetBonds()) == 2), None)
    oxygen_atom = next((atom for atom in ring_atoms if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == 0 and sum(bond.GetBondTypeAsDouble() for bond in atom.GetBonds()) == 1), None)
    if carbonyl_atom is None or oxygen_atom is None or not any(bond.GetBondTypeAsDouble() == 1 for bond in carbonyl_atom.GetBonds() if bond.GetOtherAtom(carbonyl_atom) == oxygen_atom):
        return False, "Carbonyl and oxygen not in gamma-lactone arrangement"
    
    return True, "Contains a 2-furanone skeleton with a gamma-lactone arrangement (butenolide)"