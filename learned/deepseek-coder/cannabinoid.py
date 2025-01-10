"""
Classifies: CHEBI:67194 cannabinoid
"""
"""
Classifies: CHEBI:67194 cannabinoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cannabinoid(smiles: str):
    """
    Determines if a molecule is a cannabinoid based on its SMILES string.
    Cannabinoids are characterized by either:
    1. A benzene ring with oxygen-containing functional groups and long hydrocarbon chains
    2. A fatty acid/ethanolamide structure with long hydrocarbon chains

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cannabinoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check molecular weight - cannabinoids typically have higher molecular weights
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250:
        return False, "Molecular weight too low for a cannabinoid"

    # Define patterns for cannabinoid features
    benzene_pattern = Chem.MolFromSmarts("c1ccccc1")
    oxygen_pattern = Chem.MolFromSmarts("[OX2]")
    long_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    ethanolamide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])CCO")
    fatty_acid_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2][CX4]")

    # Check for classic cannabinoid structure (benzene + oxygen + long chain)
    has_benzene = mol.HasSubstructMatch(benzene_pattern)
    has_oxygen = len(mol.GetSubstructMatches(oxygen_pattern)) >= 1
    has_long_chain = mol.HasSubstructMatch(long_chain_pattern)

    # Check for ethanolamide/fatty acid structure
    has_ethanolamide = mol.HasSubstructMatch(ethanolamide_pattern)
    has_fatty_acid = mol.HasSubstructMatch(fatty_acid_pattern)

    # Classification logic
    if (has_benzene and has_oxygen and has_long_chain):
        return True, "Contains a benzene ring with oxygen-containing functional groups and long hydrocarbon chains"
    elif (has_ethanolamide or has_fatty_acid) and has_long_chain:
        return True, "Contains ethanolamide/fatty acid structure with long hydrocarbon chain"
    elif has_benzene and has_oxygen:
        return False, "Contains benzene and oxygen but lacks characteristic long chain"
    else:
        return False, "Does not match cannabinoid structural patterns"