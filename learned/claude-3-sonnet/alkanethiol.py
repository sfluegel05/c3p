"""
Classifies: CHEBI:47908 alkanethiol
"""
"""
Classifies: alkanethiol
Definition: A compound in which a sulfanyl group, -SH, is attached to an alkyl group
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_alkanethiol(smiles: str):
    """
    Determines if a molecule is an alkanethiol based on its SMILES string.
    An alkanethiol has a sulfanyl group (-SH) attached to an alkyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkanethiol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check molecular weight - exclude large molecules like peptides
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 300:
        return False, "Molecular weight too high for simple alkanethiol"

    # Look for -SH (sulfanyl) group pattern
    sulfanyl_pattern = Chem.MolFromSmarts("[SH1]")
    if not mol.HasSubstructMatch(sulfanyl_pattern):
        return False, "No sulfanyl (-SH) group found"

    # Exclude molecules with peptide bonds
    peptide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[CX4]")
    if mol.HasSubstructMatch(peptide_pattern):
        return False, "Contains peptide bonds"

    # Count nitrogen atoms - exclude amino acids and peptides
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count > 1:
        return False, "Too many nitrogen atoms for simple alkanethiol"

    # Get all sulfur atoms
    sulfur_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16]
    
    # Check if any sulfanyl group is attached to appropriate carbon
    for s_atom in sulfur_atoms:
        neighbors = s_atom.GetNeighbors()
        if len(neighbors) == 1:  # Terminal sulfur (as in -SH)
            neighbor = neighbors[0]
            if neighbor.GetAtomicNum() == 6:  # Carbon
                # Check if carbon is sp3 or sp2 (but not aromatic)
                if (not neighbor.GetIsAromatic() and 
                    neighbor.GetHybridization() in 
                    [Chem.HybridizationType.SP3, Chem.HybridizationType.SP2]):
                    
                    # Additional check for carboxyl groups
                    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
                    if mol.HasSubstructMatch(carboxyl_pattern):
                        return False, "Contains carboxylic acid group"
                    
                    # Check for amino groups
                    amino_pattern = Chem.MolFromSmarts("[NX3;H2,H1;!$(NC=O)]")
                    if mol.HasSubstructMatch(amino_pattern):
                        return False, "Contains primary or secondary amine"
                        
                    return True, "Contains sulfanyl group (-SH) attached to an alkyl group"

    return False, "Sulfanyl group not properly attached to carbon"