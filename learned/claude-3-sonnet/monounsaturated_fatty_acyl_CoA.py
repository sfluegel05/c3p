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

    # Check for CoA moiety
    # Look for adenine pattern
    adenine_pattern = Chem.MolFromSmarts("c1nc(N)c2ncnc2n1")
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "No CoA moiety found (missing adenine)"
    
    # Look for pantetheine pattern (includes thioester)
    pantetheine_pattern = Chem.MolFromSmarts("NC(=O)CCNC(=O)CCSC(=O)")
    if not mol.HasSubstructMatch(pantetheine_pattern):
        return False, "No CoA moiety found (missing pantetheine)"

    # Count number of double bonds
    # First, get all double bonds
    double_bonds = mol.GetSubstructMatches(Chem.MolFromSmarts("C=C"))
    
    # Filter out those in adenine ring
    non_aromatic_double_bonds = []
    for bond in double_bonds:
        atom1, atom2 = mol.GetAtomWithIdx(bond[0]), mol.GetAtomWithIdx(bond[1])
        if not atom1.IsInRing() or not atom2.IsInRing():
            non_aromatic_double_bonds.append(bond)
    
    if len(non_aromatic_double_bonds) != 1:
        return False, f"Found {len(non_aromatic_double_bonds)} C=C double bonds, need exactly 1"

    # Verify thioester linkage
    thioester_pattern = Chem.MolFromSmarts("SC(=O)")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"

    # Additional check to ensure the double bond is in the fatty acyl chain
    # by checking if it's between the thioester and the end of the chain
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if thioester_matches:
        sulfur_idx = thioester_matches[0][0]
        double_bond_atoms = non_aromatic_double_bonds[0]
        
        # Create a substructure that represents the fatty acyl chain
        # by breaking at the sulfur atom
        fragments = Chem.FragmentOnBonds(mol, [mol.GetBondBetweenAtoms(sulfur_idx, 
                                      list(mol.GetAtomWithIdx(sulfur_idx).GetNeighbors())[0].GetIdx()).GetIdx()])
        
        # If fragments is None or empty, something went wrong
        if not fragments:
            return False, "Could not verify double bond position"
            
        # Check if double bond atoms are in the same fragment as thioester
        acyl_fragment = Chem.MolFromSmiles(Chem.MolToSmiles(fragments))
        if not acyl_fragment.HasSubstructMatch(Chem.MolFromSmarts("C=C")):
            return False, "Double bond not in fatty acyl chain"

    return True, "Contains CoA moiety and single double bond in fatty acyl chain"