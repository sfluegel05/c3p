"""
Classifies: CHEBI:61910 very long-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_very_long_chain_fatty_acyl_CoA(smiles: str) -> tuple[bool, str]:
    """
    Determines if a molecule is a very long-chain fatty acyl-CoA based on its SMILES string.
    Very long-chain fatty acyl-CoAs have fatty acid chains longer than C22.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple[bool, str]: (is_vlcfa_coa, reason)
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for CoA moiety
    # Look for adenine
    adenine_pattern = Chem.MolFromSmarts("c1nc(N)c2ncnc2n1")
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "No CoA moiety found (missing adenine)"
    
    # Look for thioester group (C(=O)S)
    thioester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[SX2]")
    if not mol.HasSubstructMatches(thioester_pattern):
        return False, "No thioester group found"

    # Look for phosphate groups
    phosphate_pattern = Chem.MolFromSmarts("[OX2]P(=O)([OX2])[OX2]")
    if len(mol.GetSubstructMatches(phosphate_pattern)) < 3:
        return False, "Missing phosphate groups characteristic of CoA"

    # Count carbons in the fatty acid chain
    # First, find the thioester carbon
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Could not identify thioester group"
    
    # Create a copy of the molecule to work with
    mol_copy = Chem.RWMol(mol)
    
    # Break the molecule at the thioester sulfur to isolate the fatty acid portion
    thioester_carbon = thioester_matches[0][0]
    bonds_to_break = []
    for bond in mol_copy.GetBonds():
        if bond.GetBeginAtomIdx() == thioester_carbon and bond.GetEndAtom().GetAtomicNum() == 16:  # Sulfur
            bonds_to_break.append(bond.GetIdx())
        elif bond.GetEndAtomIdx() == thioester_carbon and bond.GetBeginAtom().GetAtomicNum() == 16:  # Sulfur
            bonds_to_break.append(bond.GetIdx())
    
    # Break bonds in reverse order to maintain indices
    for bond_idx in sorted(bonds_to_break, reverse=True):
        mol_copy.RemoveBond(mol_copy.GetBondWithIdx(bond_idx).GetBeginAtomIdx(),
                           mol_copy.GetBondWithIdx(bond_idx).GetEndAtomIdx())
    
    # Get the fragment containing the fatty acid
    fragments = Chem.GetMolFrags(mol_copy, asMols=True)
    fatty_acid = None
    for frag in fragments:
        if frag.HasSubstructMatch(Chem.MolFromSmarts("[CX3](=[OX1])")):
            fatty_acid = frag
            break
    
    if fatty_acid is None:
        return False, "Could not isolate fatty acid portion"
    
    # Count carbons in fatty acid portion
    carbon_count = sum(1 for atom in fatty_acid.GetAtoms() if atom.GetAtomicNum() == 6)
    
    if carbon_count <= 22:
        return False, f"Fatty acid chain length (C{carbon_count}) not greater than C22"
    
    # Additional checks for fatty acid characteristics
    if carbon_count > 40:
        return False, f"Fatty acid chain length (C{carbon_count}) unreasonably long"

    return True, f"Very long-chain fatty acyl-CoA with C{carbon_count} fatty acid chain"