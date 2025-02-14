"""
Classifies: CHEBI:73702 wax
"""
"""
Classifies: CHEBI:82198 wax
"""

from rdkit import Chem

def is_wax(smiles: str):
    """
    Determines if a molecule is a wax based on its SMILES string.
    A wax is defined as an organic compound or mixture of compounds
    that is composed of long-chain molecules and is malleable at ambient temperatures.
    Commonly, waxes are esters formed from long-chain fatty acids and long-chain alcohols.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a wax, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify ester bonds: bond between carbonyl carbon and ester oxygen
    ester_bond_indices = []
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        if (atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 8) or (atom1.GetAtomicNum() == 8 and atom2.GetAtomicNum() == 6):
            # Check if carbon is carbonyl carbon (has double bond to oxygen)
            if atom1.GetAtomicNum() == 6:
                carbon = atom1
                oxygen = atom2
            else:
                carbon = atom2
                oxygen = atom1
            is_carbonyl = False
            for nbr in carbon.GetNeighbors():
                if nbr.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(carbon.GetIdx(), nbr.GetIdx()).GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    is_carbonyl = True
                    break
            # Check if oxygen is connected to another carbon
            if is_carbonyl:
                for nbr in oxygen.GetNeighbors():
                    if nbr.GetIdx() != carbon.GetIdx() and nbr.GetAtomicNum() == 6:
                        ester_bond_indices.append(bond.GetIdx())
                        break

    if not ester_bond_indices:
        return False, "No ester groups found"

    # Define threshold for long-chain (number of carbons)
    chain_length_threshold = 12

    # Break the molecule at ester bonds and analyze fragments
    frags = Chem.FragmentOnBonds(mol, ester_bond_indices, addDummies=False)
    frag_mols = Chem.GetMolFrags(frags, asMols=True, sanitizeFrags=True)

    # Count carbons in each fragment
    long_chain_fragments = []
    for frag in frag_mols:
        c_count = sum(1 for atom in frag.GetAtoms() if atom.GetAtomicNum() == 6)
        if c_count >= chain_length_threshold:
            long_chain_fragments.append(frag)

    # Check if there are at least two long-chain fragments (from ester bond)
    if len(long_chain_fragments) >= 2:
        return True, f"Ester with long chains found (number of long chains: {len(long_chain_fragments)})"
    else:
        return False, "No ester with long chains found"