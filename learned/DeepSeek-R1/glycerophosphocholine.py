"""
Classifies: CHEBI:36313 glycerophosphocholine
"""
"""
Classifies: CHEBI:17585 glycerophosphocholine
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_glycerophosphocholine(smiles: str):
    """
    Determines if a molecule is a glycerophosphocholine based on its SMILES string.
    A glycerophosphocholine consists of a glycerol backbone with a phosphocholine group and two fatty acid chains.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycerophosphocholine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for phosphocholine group (P connected to choline and glycerol oxygen)
    phosphocholine_pattern = Chem.MolFromSmarts("[O]-P(=O)([O-])-O-C-C-[N+](C)(C)C")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "Phosphocholine group not found"

    # Find phosphorus atom in phosphocholine group
    matches = mol.GetSubstructMatches(phosphocholine_pattern)
    p_atom = matches[0][1]  # P is the second atom in the SMARTS pattern [O]-P...

    # Get oxygen connecting P to glycerol (the first O in the pattern)
    connecting_o = matches[0][0]
    glycerol_carbon = [x.GetIdx() for x in mol.GetAtomWithIdx(connecting_o).GetNeighbors() if x.GetIdx() != p_atom][0]

    # Check that glycerol carbon has exactly two other oxygen attachments (for fatty acid chains)
    glycerol_c = mol.GetAtomWithIdx(glycerol_carbon)
    o_neighbors = [n for n in glycerol_c.GetNeighbors() if n.GetAtomicNum() == 8 and n.GetIdx() != connecting_o]
    if len(o_neighbors) != 2:
        return False, f"Glycerol carbon has {len(o_neighbors)} oxygen attachments (need 2)"

    # Check each oxygen is part of ester or ether group
    valid_chain_count = 0
    for o in o_neighbors:
        # Check for ester (O-C=O) or ether (O-C without carbonyl)
        ester_match = False
        ether_match = False
        
        # Check ester: O connected to C=O
        for neighbor in o.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:
                for bond in neighbor.GetBonds():
                    if bond.GetBondType() == Chem.BondType.DOUBLE and bond.GetBeginAtomIdx() == neighbor.GetIdx() and bond.GetEndAtom().GetAtomicNum() == 8:
                        ester_match = True
                        break
        # Check ether: O connected to C (non-carbonyl)
        if not ester_match:
            for neighbor in o.GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and any(b.GetBondType() != Chem.BondType.DOUBLE for b in neighbor.GetBonds() if b.GetOtherAtomIdx(neighbor.GetIdx()) != o.GetIdx()):
                    ether_match = True
                    break
        
        if ester_match or ether_match:
            valid_chain_count += 1

    if valid_chain_count < 2:
        return False, f"Found {valid_chain_count} valid fatty acid chains (need 2)"

    # Optional: Check molecular weight to filter very small molecules
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt < 400:
        return False, f"Molecular weight {mol_wt:.1f} too low for typical glycerophosphocholine"

    return True, "Glycerol backbone with phosphocholine and two fatty acid chains"