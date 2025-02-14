"""
Classifies: CHEBI:63534 monoamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monoamine(smiles: str):
    """
    Determines if a molecule is a monoamine based on its SMILES string.
    A monoamine is defined as an arylalkylamino compound which contains one amino group connected to an aromatic ring by a two-carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoamine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for exactly one amino group (primary, secondary, tertiary or quarternary, including protonated forms)
    amino_pattern = Chem.MolFromSmarts("[NX3;H2,H1,H0;!$(N=*)]") # Primary, secondary or tertiary amino groups
    ammonium_pattern = Chem.MolFromSmarts("[NX4;!$(N=*)]") # Quarternary ammonium
    amino_matches = mol.GetSubstructMatches(amino_pattern)
    ammonium_matches = mol.GetSubstructMatches(ammonium_pattern)
    
    if len(amino_matches) + len(ammonium_matches) != 1:
        return False, f"Molecule must have exactly one amino group. Found {len(amino_matches) + len(ammonium_matches)}"
    
    # Check for the presence of at least one substituted benzene ring
    aromatic_ring_pattern = Chem.MolFromSmarts("c1ccccc1")
    if not mol.HasSubstructMatch(aromatic_ring_pattern):
        return False, "Molecule must contain at least one aromatic ring"
    
    # Check for two carbon chain linking amino group to aromatic ring
    two_carbon_chain_pattern = Chem.MolFromSmarts("[NX3,NX4][CX4][CX4]c1ccccc1")
    if not mol.HasSubstructMatch(two_carbon_chain_pattern):
       
        # Check for variants with a non-aromatic chain (not the preferred case but it can happen)
        two_carbon_chain_pattern_variant = Chem.MolFromSmarts("[NX3,NX4][CX4][CX4][CX3,CX4]") # Two carbons connecting to a non-aromatic carbon
        if not mol.HasSubstructMatch(two_carbon_chain_pattern_variant):
            return False, "Amino group must be connected to the aromatic ring by a two-carbon chain."

    
    return True, "Molecule contains one amino group attached to an aromatic ring via a two-carbon chain"