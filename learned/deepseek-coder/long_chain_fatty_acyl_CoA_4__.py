"""
Classifies: CHEBI:83139 long-chain fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA(4-) based on its SMILES string.
    A long-chain fatty acyl-CoA(4-) is a molecule with a long-chain fatty acid (14-24 carbons)
    attached to a CoA moiety via a thioester bond, with deprotonated phosphate and diphosphate groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for CoA moiety pattern (deprotonated phosphate and diphosphate groups)
    coa_pattern = Chem.MolFromSmarts("[O-]P(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA moiety with deprotonated phosphate and diphosphate groups found"

    # Check for thioester bond (S-C(=O)-) connected to the CoA moiety
    thioester_pattern = Chem.MolFromSmarts("[SX2][CX3](=[OX1])")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if len(thioester_matches) == 0:
        return False, "No thioester bond found"

    # Ensure the thioester bond is connected to the CoA moiety
    coa_atoms = set(mol.GetSubstructMatch(coa_pattern))
    thioester_connected = False
    for match in thioester_matches:
        if any(atom in coa_atoms for atom in match):
            thioester_connected = True
            break
    if not thioester_connected:
        return False, "Thioester bond not connected to CoA moiety"

    # Count carbons in the fatty acid chain
    # Find the carbon chain connected to the thioester bond
    fatty_acid_carbons = set()
    for match in thioester_matches:
        if any(atom in coa_atoms for atom in match):
            # The carbon in the thioester bond is part of the fatty acid chain
            fatty_acid_carbons.add(match[1])
            # Traverse the chain from this carbon
            stack = [match[1]]
            while stack:
                atom_idx = stack.pop()
                atom = mol.GetAtomWithIdx(atom_idx)
                if atom.GetAtomicNum() == 6:
                    fatty_acid_carbons.add(atom_idx)
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetIdx() not in fatty_acid_carbons and neighbor.GetAtomicNum() == 6:
                            stack.append(neighbor.GetIdx())

    # Check if the fatty acid chain length is within the expected range (14-24 carbons)
    if len(fatty_acid_carbons) < 14 or len(fatty_acid_carbons) > 24:
        return False, f"Fatty acid chain length ({len(fatty_acid_carbons)} carbons) is not within the expected range (14-24 carbons)"

    # Check molecular weight - long-chain fatty acyl-CoA(4-) typically >700 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 700:
        return False, "Molecular weight too low for long-chain fatty acyl-CoA(4-)"

    return True, "Contains a long-chain fatty acid attached to a CoA moiety via a thioester bond, with deprotonated phosphate and diphosphate groups"