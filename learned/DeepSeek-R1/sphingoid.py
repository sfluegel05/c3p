"""
Classifies: CHEBI:35785 sphingoid
"""
"""
Classifies: CHEBI: ??? sphingoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sphingoid(smiles: str):
    """
    Determines if a molecule is a sphingoid based on its SMILES string.
    Sphingoids are characterized by a long-chain amino alcohol backbone with hydroxyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sphingoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for amino group (NH2, NH3+, or amide)
    amino_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:
            # Check for primary/secondary amine or ammonium
            if atom.GetTotalDegree() <= 3 and (atom.GetFormalCharge() in (0, 1) or atom.GetTotalNumHs() >= 1):
                amino_found = True
                break
            # Check for amide (N connected to carbonyl)
            if any(bond.GetBondType() == Chem.BondType.SINGLE 
                   and bond.GetOtherAtom(atom).GetAtomicNum() == 6 
                   and any(nbr.GetAtomicNum() == 8 for nbr in bond.GetOtherAtom(atom).GetNeighbors())
                   for bond in atom.GetBonds()):
                amino_found = True
                break

    if not amino_found:
        return False, "No amino group detected"

    # Check for at least one hydroxyl group adjacent to amino-bearing carbon
    hydroxyl_adjacent = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:
                    for bond in neighbor.GetBonds():
                        other = bond.GetOtherAtom(neighbor)
                        if other.GetAtomicNum() == 8 and other.GetTotalNumHs() >= 1:
                            hydroxyl_adjacent = True
                            break
                    if hydroxyl_adjacent:
                        break
            if hydroxyl_adjacent:
                break

    if not hydroxyl_adjacent:
        return False, "No hydroxyl group adjacent to amino group"

    # Check for long aliphatic chain (>=12 carbons)
    chain_length = rdMolDescriptors.CalcNumAcyclicStereocenters(mol) + 1  # Approximation
    if chain_length < 12:
        # Fallback to exact chain length calculation
        chain_length = max(
            len(Chem.GetLongestChain(mol, includeDoubleBonds=True)),
            len(Chem.GetLongestChain(mol, includeDoubleBonds=False))
        )
    
    if chain_length < 12:
        return False, f"Chain length {chain_length} < 12"

    return True, "Long-chain amino alcohol with adjacent hydroxyl group"