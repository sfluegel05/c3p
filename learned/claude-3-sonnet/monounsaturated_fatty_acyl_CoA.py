"""
Classifies: CHEBI:139575 monounsaturated fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monounsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acyl-CoA based on its SMILES string.
    Must have exactly one carbon-carbon double bond in the fatty acyl chain.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) - (is_monounsaturated_fatty_acyl_CoA, reason)
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for adenine nucleotide (part of CoA)
    adenine_pattern = Chem.MolFromSmarts("c1nc(N)c2ncnc2n1")
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "No CoA moiety found (missing adenine)"
    
    # Look for phosphate groups characteristic of CoA
    phosphate_pattern = Chem.MolFromSmarts("OP(=O)(O)OP(=O)(O)O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "Missing characteristic CoA phosphate groups"

    # Look for thioester linkage
    thioester_pattern = Chem.MolFromSmarts("[C;X3](=[O;X1])[S;X2]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester linkage found"

    # Count all double bonds in the molecule
    double_bond_pattern = Chem.MolFromSmarts("[C;X3]=[C;X3]")
    all_double_bonds = mol.GetSubstructMatches(double_bond_pattern)
    
    # Filter out double bonds in adenine ring
    non_aromatic_double_bonds = []
    for bond in all_double_bonds:
        atom1 = mol.GetAtomWithIdx(bond[0])
        atom2 = mol.GetAtomWithIdx(bond[1])
        if not (atom1.GetIsAromatic() or atom2.GetIsAromatic()):
            non_aromatic_double_bonds.append(bond)

    if len(non_aromatic_double_bonds) != 1:
        return False, f"Found {len(non_aromatic_double_bonds)} non-aromatic C=C double bonds, need exactly 1"

    # Verify the double bond is in a fatty acid chain by checking connectivity to thioester
    thioester_sulfur = mol.GetAtomWithIdx(thioester_matches[0][2])  # Get sulfur atom
    double_bond_atoms = [mol.GetAtomWithIdx(non_aromatic_double_bonds[0][0]),
                        mol.GetAtomWithIdx(non_aromatic_double_bonds[0][1])]
    
    # Check if there's a path between thioester and double bond
    path_exists = False
    for atom in double_bond_atoms:
        if len(Chem.GetShortestPath(mol, thioester_sulfur.GetIdx(), atom.GetIdx())) > 0:
            path_exists = True
            break
    
    if not path_exists:
        return False, "Double bond not connected to fatty acyl chain"

    # Count carbons in the fatty acid chain
    carbon_chain = Chem.MolFromSmarts("[C;X4,X3]~[C;X4,X3]~[C;X4,X3]~[C;X4,X3]")
    if not mol.HasSubstructMatch(carbon_chain):
        return False, "Fatty acid chain too short"

    # Additional check for pantetheine arm
    pantetheine = Chem.MolFromSmarts("NCCC(=O)NCCS")
    if not mol.HasSubstructMatch(pantetheine):
        return False, "Missing pantetheine arm of CoA"

    return True, "Contains CoA moiety with single double bond in fatty acyl chain"