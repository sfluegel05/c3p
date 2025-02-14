"""
Classifies: CHEBI:35785 sphingoid
"""
"""
Classifies: CHEBI:18089 sphingoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sphingoid(smiles: str):
    """
    Determines if a molecule is a sphingoid based on its SMILES string.
    A sphingoid is defined as 'Sphinganine, its homologs and stereoisomers, and the hydroxy and unsaturated derivatives of these compounds'.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sphingoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for sphingoid backbone pattern (long alkyl chain with terminal amino alcohol)
    backbone_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[NX3][CX4][CX3][OX2H]")
    if not mol.HasSubstructMatch(backbone_pattern):
        return False, "No sphingoid backbone found"

    # Look for common modifications
    glucosyl_pattern = Chem.MolFromSmarts("[OX2]C[C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O")
    phosphocholine_pattern = Chem.MolFromSmarts("[NX4+]([C])(C)[C]OP(=O)([O-])[O-]")

    # Check for unsaturated and branched alkyl chains
    unsaturated_chain_pattern = Chem.MolFromSmarts("[CX3]=&!@[CX3]")
    branched_chain_pattern = Chem.MolFromSmarts("[CX4]([CX4])([CX4])[CX4]")

    # Check for common properties
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)

    # Build reason for classification
    reason = "Contains sphingoid backbone: "
    if mol.HasSubstructMatch(glucosyl_pattern):
        reason += "with glucosyl group, "
    if mol.HasSubstructMatch(phosphocholine_pattern):
        reason += "with phosphocholine group, "
    if mol.HasSubstructMatch(unsaturated_chain_pattern):
        reason += "with unsaturated alkyl chain, "
    if mol.HasSubstructMatch(branched_chain_pattern):
        reason += "with branched alkyl chain, "
    reason += f"mol. wt. {mol_wt:.2f} Da, {c_count} carbons, {o_count} oxygens, {n_count} nitrogens, {n_rotatable} rotatable bonds"

    # Check additional criteria
    if mol_wt < 250 or mol_wt > 800:
        return False, "Molecular weight outside typical range for sphingoids"
    if c_count < 12 or c_count > 30:
        return False, "Carbon count outside typical range for sphingoids"
    if o_count < 2 or o_count > 8:
        return False, "Oxygen count outside typical range for sphingoids"
    if n_count != 1:
        return False, "Must have exactly 1 nitrogen"
    if n_rotatable < 10:
        return False, "Not enough rotatable bonds for long alkyl chain"

    return True, reason