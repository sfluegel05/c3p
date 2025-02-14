"""
Classifies: CHEBI:90546 medium-chain fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_medium_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA(4-) based on its SMILES string.
    A medium-chain fatty acyl-CoA(4-) is an acyl-CoA molecule with an acyl chain of 6-12 carbons,
    in which the phosphate and diphosphate groups are deprotonated (negative charges).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a medium-chain fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for phosphorus atoms (CoA derivative)
    contains_phosphorus = any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms())
    if not contains_phosphorus:
        return False, "No phosphorus atoms found; molecule is not a CoA derivative"

    # Check for adenine moiety (part of CoA)
    adenine_pattern = Chem.MolFromSmarts('n1cnc2c1ncnc2N')
    adenine_matches = mol.GetSubstructMatches(adenine_pattern)
    if not adenine_matches:
        return False, "No adenine moiety found; molecule is not a CoA derivative"

    # Look for thioester linkage C(=O)-S
    thioester_pattern = Chem.MolFromSmarts('C(=O)S')
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester linkage found"

    # Use the first thioester match
    acyl_carbon_idx = thioester_matches[0][0]  # Carbonyl carbon
    sulfur_idx = thioester_matches[0][2]       # Sulfur atom

    # Break bond between acyl carbon and sulfur to isolate acyl chain
    bond = mol.GetBondBetweenAtoms(acyl_carbon_idx, sulfur_idx)
    if bond is None:
        return False, "Cannot find bond between acyl carbon and sulfur"

    bond_idx = bond.GetIdx()
    mol_frag = Chem.FragmentOnBonds(mol, [bond_idx])

    # Get fragments
    frags = Chem.GetMolFrags(mol_frag, asMols=True, sanitizeFrags=True)
    acyl_chain_frag = None
    for frag in frags:
        has_sulfur = any(atom.GetAtomicNum() == 16 for atom in frag.GetAtoms())
        if not has_sulfur:
            acyl_chain_frag = frag
            break

    if acyl_chain_frag is None:
        return False, "Could not isolate acyl chain fragment"

    # Count carbons in acyl chain
    carbon_count = sum(1 for atom in acyl_chain_frag.GetAtoms() if atom.GetAtomicNum() == 6)

    if 6 <= carbon_count <= 12:
        return True, f"Acyl chain has {carbon_count} carbons, which is within medium-chain length (6-12)"
    else:
        return False, f"Acyl chain has {carbon_count} carbons, not within medium-chain length (6-12)"