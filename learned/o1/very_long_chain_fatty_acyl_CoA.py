"""
Classifies: CHEBI:61910 very long-chain fatty acyl-CoA
"""
"""
Classifies: very long-chain fatty acyl-CoA
Definition: A fatty acyl-CoA in which the fatty acyl group has a chain length greater than C22.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_very_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acyl-CoA based on its SMILES string.
    A very long-chain fatty acyl-CoA is a fatty acyl-CoA molecule where the acyl chain has more than 22 carbons.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a very long-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find the thioester bond: C(=O)-S
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    matches = mol.GetSubstructMatches(thioester_pattern)
    if not matches:
        return False, "No thioester linkage found"

    # Assume the first match is the acyl-CoA linkage
    match = matches[0]
    carbonyl_c_idx = match[0]
    sulfur_idx = match[2]  # Index of the sulfur atom

    # Get the bond index of the C-S bond to break
    cs_bond = mol.GetBondBetweenAtoms(carbonyl_c_idx, sulfur_idx)
    if cs_bond is None:
        return False, "Failed to identify thioester bond"

    bond_idx = cs_bond.GetIdx()

    # Break the thioester bond to separate acyl chain from CoA
    mol_frag = Chem.FragmentOnBonds(mol, [bond_idx])
    frags = Chem.GetMolFrags(mol_frag, asMols=True, sanitizeFrags=True)

    # Identify the acyl chain fragment
    acyl_chain = None
    for frag in frags:
        # Check if fragment contains the carbonyl carbon (part of acyl chain)
        if frag.HasSubstructMatch(Chem.MolFromSmarts("C(=O)[CH0,CH1,CH2,CH3]")):
            acyl_chain = frag
            break

    if acyl_chain is None:
        return False, "Failed to isolate acyl chain"

    # Count the number of carbon atoms in the acyl chain
    num_carbons = sum(1 for atom in acyl_chain.GetAtoms() if atom.GetAtomicNum() == 6)

    if num_carbons > 22:
        return True, f"Acyl chain has {num_carbons} carbons, which is greater than 22"
    else:
        return False, f"Acyl chain has {num_carbons} carbons, which is not greater than 22"