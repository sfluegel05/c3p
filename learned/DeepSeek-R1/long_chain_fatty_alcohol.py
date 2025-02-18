"""
Classifies: CHEBI:17135 long-chain fatty alcohol
"""
"""
Classifies: long-chain fatty alcohol (C13-C22)
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a long-chain fatty alcohol (C13-C22 chain) based on its SMILES string.
    A fatty alcohol has a hydroxyl group attached to an aliphatic carbon chain of 13-22 carbons.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Must have at least one hydroxyl group
    if not any(atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() == 1 for atom in mol.GetAtoms()):
        return False, "No hydroxyl group"

    # Create carbon-only structure with single bonds for chain analysis
    emol = Chem.RWMol()
    carbon_map = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon only
            new_atom = Chem.Atom(6)
            idx = emol.AddAtom(new_atom)
            carbon_map[atom.GetIdx()] = idx

    # Add single bonds between carbons
    for bond in mol.GetBonds():
        if bond.GetBeginAtom().GetAtomicNum() == 6 and bond.GetEndAtom().GetAtomicNum() == 6:
            i = carbon_map[bond.GetBeginAtom().GetIdx()]
            j = carbon_map[bond.GetEndAtom().GetIdx()]
            emol.AddBond(i, j, Chem.BondType.SINGLE)

    # Calculate longest chain in carbon skeleton
    mol_carbons = emol.GetMol()
    if mol_carbons.GetNumAtoms() == 0:
        return False, "No carbon chain"
    
    longest_chain = rdMolDescriptors.CalcLongestChain(mol_carbons)

    # Verify chain length requirement
    if 13 <= longest_chain <= 22:
        return True, f"Contains C{longest_chain} chain with hydroxyl group"
    return False, f"Chain length {longest_chain} not in 13-22 range"