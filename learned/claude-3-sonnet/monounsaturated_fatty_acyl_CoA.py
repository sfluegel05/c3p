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

    # Neutralize charges for consistent handling
    uncharger = Chem.AllChem.UnchargeModel()
    uncharged_mol = uncharger.uncharge(mol)

    # Check for complete CoA structure
    # Look for adenine nucleotide
    adenine_pattern = Chem.MolFromSmarts("c1nc(N)c2ncnc2n1")
    if not uncharged_mol.HasSubstructMatch(adenine_pattern):
        return False, "No CoA moiety found (missing adenine)"
    
    # Look for ribose-phosphate backbone
    ribose_phosphate = Chem.MolFromSmarts("[CH2]OP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(=O)(O)O)n1cnc2c(N)ncnc12")
    if not uncharged_mol.HasSubstructMatch(ribose_phosphate):
        return False, "Missing or incorrect ribose-phosphate backbone"

    # Look for pantetheine arm with correct structure
    pantetheine = Chem.MolFromSmarts("CC(C)(COP(=O)(O))[C@@H](O)C(=O)NCCC(=O)NCCSC(=O)")
    if not uncharged_mol.HasSubstructMatch(pantetheine):
        return False, "Missing or incorrect pantetheine arm"

    # Look for thioester linkage specifically connected to a carbon chain
    thioester_pattern = Chem.MolFromSmarts("[S;X2][C;X3](=[O;X1])[C;X4]")
    thioester_matches = uncharged_mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester linkage found or incorrect connectivity"

    # Count non-aromatic, non-conjugated double bonds in the fatty acid portion
    double_bond_pattern = Chem.MolFromSmarts("[CX4,CX3][CX3]=[CX3][CX4,CX3]")
    double_bonds = uncharged_mol.GetSubstructMatches(double_bond_pattern)
    
    # Filter out double bonds in rings or conjugated systems
    valid_double_bonds = []
    for bond in double_bonds:
        atoms = [uncharged_mol.GetAtomWithIdx(i) for i in bond]
        # Check if any atom is in a ring
        if any(atom.IsInRing() for atom in atoms):
            continue
        # Check for conjugation
        if any(atom.GetIsConjugated() for atom in atoms):
            continue
        valid_double_bonds.append(bond)

    if len(valid_double_bonds) != 1:
        return False, f"Found {len(valid_double_bonds)} valid C=C double bonds, need exactly 1"

    # Verify fatty acid chain length and connectivity
    thioester_sulfur = uncharged_mol.GetAtomWithIdx(thioester_matches[0][0])
    carbon_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    carbon_chains = uncharged_mol.GetSubstructMatches(carbon_chain_pattern)
    
    # Check if there's a continuous carbon chain from thioester to double bond
    valid_chain = False
    for chain in carbon_chains:
        if any(idx in chain for db_bond in valid_double_bonds for idx in db_bond):
            path_to_thioester = False
            for atom_idx in chain:
                if len(Chem.FindAllPathsOfLengthN(uncharged_mol, 3, thioester_sulfur.GetIdx(), atom_idx)) > 0:
                    path_to_thioester = True
                    break
            if path_to_thioester:
                valid_chain = True
                break

    if not valid_chain:
        return False, "No valid fatty acid chain containing the double bond"

    # Additional check for minimum chain length from thioester
    chain_length = max(len(path) for path in Chem.FindAllPathsOfLengthN(uncharged_mol, 8, thioester_sulfur.GetIdx(), -1))
    if chain_length < 6:
        return False, "Fatty acid chain too short"

    return True, "Contains CoA moiety and single double bond in fatty acyl chain"