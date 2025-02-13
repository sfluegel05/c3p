"""
Classifies: CHEBI:139575 monounsaturated fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monounsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acyl-CoA based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) - (is_monounsaturated_fatty_acyl_CoA, reason)
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for CoA moiety components
    # Look for adenine pattern
    adenine_pattern = Chem.MolFromSmarts("c1nc(N)c2ncnc2n1")
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "No CoA moiety found (missing adenine)"
    
    # Look for phosphate groups (need at least 3 for CoA)
    phosphate_pattern = Chem.MolFromSmarts("OP(=O)(O)O")
    phosphate_matches = len(mol.GetSubstructMatches(phosphate_pattern))
    if phosphate_matches < 3:
        return False, f"Insufficient phosphate groups for CoA (found {phosphate_matches}, need ≥3)"

    # Look for thioester linkage (more flexible pattern)
    thioester_pattern = Chem.MolFromSmarts("[S;X2][C;X3](=[O;X1])")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"

    # Look for amide bonds characteristic of pantetheine
    amide_pattern = Chem.MolFromSmarts("[NX3][C;X3](=[O;X1])[CX4]")
    amide_matches = len(mol.GetSubstructMatches(amide_pattern))
    if amide_matches < 2:
        return False, f"Missing pantetheine amide bonds (found {amide_matches}, need ≥2)"

    # Count non-aromatic double bonds
    double_bond_pattern = Chem.MolFromSmarts("[CX3]=[CX3]")
    double_bonds = mol.GetSubstructMatches(double_bond_pattern)
    
    # Filter out aromatic bonds
    non_aromatic_double_bonds = []
    for bond in double_bonds:
        atom1, atom2 = mol.GetAtomWithIdx(bond[0]), mol.GetAtomWithIdx(bond[1])
        if not (atom1.IsInRing() and atom2.IsInRing() and atom1.GetIsAromatic() and atom2.GetIsAromatic()):
            non_aromatic_double_bonds.append(bond)
    
    if len(non_aromatic_double_bonds) != 1:
        return False, f"Found {len(non_aromatic_double_bonds)} C=C double bonds, need exactly 1"

    # Verify the presence of a fatty acid chain
    # Get carbon chain length
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 20:  # CoA itself has about 23 carbons, so fatty acyl-CoA should have more
        return False, f"Carbon count ({carbon_count}) too low for fatty acyl-CoA"

    # Verify the double bond is in the fatty acyl portion
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if thioester_matches:
        sulfur_idx = thioester_matches[0][0]
        for db_atoms in non_aromatic_double_bonds:
            # Check if we can find a path from the sulfur to both double bond carbons
            # that doesn't go through the CoA part
            paths_exist = any(len(Chem.FindAllPathsOfLengthN(mol, 2, sulfur_idx, db_atom)) > 0 
                            for db_atom in db_atoms)
            if not paths_exist:
                return False, "Double bond not in fatty acyl chain"

    return True, "Contains CoA moiety and single double bond in fatty acyl chain"