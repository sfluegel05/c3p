"""
Classifies: CHEBI:17135 long-chain fatty alcohol
"""
"""
Classifies: CHEBI:26766 long-chain fatty alcohol
A fatty alcohol with a chain length ranging from C13 to C22.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a long-chain fatty alcohol based on its SMILES string.
    A long-chain fatty alcohol has a carbon chain length ranging from C13 to C22,
    with an -OH group attached directly or via a short linker.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for -OH group
    alcohol_pattern = Chem.MolFromSmarts("[OX1H]")
    alcohol_matches = mol.GetSubstructMatches(alcohol_pattern)
    if not alcohol_matches:
        return False, "No alcohol (-OH) group found"

    # Check for long carbon chain (C13-C22)
    chain_pattern = Chem.MolFromSmarts("[CX4]([CX4])[CX4]([CX4])[CX4]([CX4])[CX4]([CX4])[CX4]([CX4])[CX4]([CX4])[CX4]([CX4])[CX4]([CX4])[CX4]([CX4])[CX4]([CX4])[CX4]([CX4])[CX4]([CX4])[CX4]([CX4])")
    chain_matches = mol.GetSubstructMatches(chain_pattern)
    if not chain_matches:
        return False, "No long carbon chain (C13-C22) found"

    # Check for linear, non-cyclic structure
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains cyclic structures"

    # Check for connectivity between alcohol group and carbon chain
    for alcohol_idx in alcohol_matches:
        alcohol_atom = mol.GetAtomWithIdx(alcohol_idx)
        for chain_idx in chain_matches:
            chain_atom = mol.GetAtomWithIdx(chain_idx)
            if mol.GetBondBetweenAtoms(alcohol_atom.GetIdx(), chain_atom.GetIdx()) is not None:
                break
        else:
            continue
        break
    else:
        return False, "Alcohol group not connected to carbon chain"

    # Check molecular weight range (typically 200-350 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 350:
        return False, "Molecular weight outside typical range for long-chain fatty alcohols"

    return True, "Meets criteria for long-chain fatty alcohol (C13-C22 chain, -OH group)"