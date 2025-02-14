"""
Classifies: CHEBI:139575 monounsaturated fatty acyl-CoA
"""
"""
Classifies: CHEBI:36691 monounsaturated fatty acyl-CoA
"""

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monounsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acyl-CoA based on its SMILES string.
    A monounsaturated fatty acyl-CoA is an unsaturated fatty acyl-CoA with one carbon-carbon double bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monounsaturated fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Load reference CoA molecule
    ref_coa_mol = Chem.MolFromSmiles("C(C)(COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(=O)([O-])=O)n1cnc2c(N)ncnc12)[C@@H](O)C(=O)NCCC(=O)NCCS")
    
    # Check for CoA moiety using molecular fingerprint similarity
    coa_fingerprint = AllChem.GetMorganFingerprint(ref_coa_mol, 2)
    mol_fingerprint = AllChem.GetMorganFingerprint(mol, 2)
    coa_similarity = DataStructs.TanimotoSimilarity(coa_fingerprint, mol_fingerprint)
    if coa_similarity < 0.8:
        return False, "No CoA moiety found"
    
    # Look for fatty acyl chain
    alkyl_pattern = Chem.MolFromSmarts("[CX4]([CX4H3])([CX4H3])([CX4H3])[CX4H2]")
    if not mol.HasSubstructMatch(alkyl_pattern):
        return False, "No fatty acyl chain found"
    
    # Count double bonds
    num_double_bonds = Chem.rdMolDescriptors.CalcNumRotatableBonds(mol) - sum(1 for b in mol.GetBonds() if b.GetBondType() == Chem.BondType.SINGLE)
    if num_double_bonds != 1:
        return False, f"Found {num_double_bonds} double bonds, expected 1"
    
    # Count carbon atoms in fatty acyl chain
    chain_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == "C" and atom.GetIsAromatic() == 0]
    chain_length = sum(1 for atom in chain_atoms if atom.GetDegree() <= 2)
    if chain_length < 4:
        return False, "Fatty acyl chain too short"
    
    return True, "Contains a CoA moiety and a monounsaturated fatty acyl chain"