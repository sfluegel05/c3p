"""
Classifies: CHEBI:26935 tetraterpenoid
"""
"""
Classifies: CHEBI:26964 tetraterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tetraterpenoid(smiles: str):
    """
    Determines if a molecule is a tetraterpenoid based on its SMILES string.
    A tetraterpenoid is a terpenoid derived from a tetraterpene (C40 skeleton),
    which may be rearranged or modified by the removal of skeletal atoms (generally methyl groups).
    It typically features a long conjugated polyene chain derived from eight isoprene units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetraterpenoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 30 or c_count > 50:
        return False, f"Carbon count ({c_count}) not in range for tetraterpenoid (30-50 carbons)"
    
    # Calculate exact molecular weight (approximate range for tetraterpenoids)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500 or mol_wt > 700:
        return False, f"Molecular weight ({mol_wt:.2f}) not in range for tetraterpenoid (500-700 Da)"
    
    # Find the longest conjugated polyene chain
    # We consider a chain of alternating single and double bonds between carbons
    mol = Chem.AddHs(mol)
    longest_chain = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            if atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 6:
                # Start traversing the conjugated system
                visited_bonds = set()
                stack = [(bond, 1)]
                while stack:
                    current_bond, length = stack.pop()
                    visited_bonds.add(current_bond.GetIdx())
                    next_atoms = [current_bond.GetBeginAtom(), current_bond.GetEndAtom()]
                    for atom in next_atoms:
                        for neighbor_bond in atom.GetBonds():
                            if neighbor_bond.GetIdx() == current_bond.GetIdx():
                                continue
                            if neighbor_bond.GetIdx() in visited_bonds:
                                continue
                            if neighbor_bond.GetBondType() in [Chem.rdchem.BondType.SINGLE, Chem.rdchem.BondType.DOUBLE]:
                                # Check for alternating bond types
                                if neighbor_bond.GetBondType() != current_bond.GetBondType():
                                    stack.append((neighbor_bond, length + 1))
                if length > longest_chain:
                    longest_chain = length
    if longest_chain < 18:
        return False, f"Longest conjugated polyene chain length is {longest_chain}, need at least 18"

    # Check for presence of multiple methyl groups (CH3)
    methyl_group = Chem.MolFromSmarts("[CH3]")
    methyl_matches = mol.GetSubstructMatches(methyl_group)
    if len(methyl_matches) < 6:
        return False, f"Found {len(methyl_matches)} methyl groups, need at least 6"
    
    # Optional: Check for isoprene units flexibility
    # Approximate number of isoprene units
    isoprene_units = c_count / 5
    if isoprene_units < 6 or isoprene_units > 10:
        return False, f"Number of isoprene units ({isoprene_units:.1f}) not in range for tetraterpenoid (6-10 units)"
    
    return True, "Molecule meets criteria for tetraterpenoid (long conjugated polyene chain derived from C40 skeleton)"