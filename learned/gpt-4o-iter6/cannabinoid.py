"""
Classifies: CHEBI:67194 cannabinoid
"""
"""
Classifies: cannabinoids
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cannabinoid(smiles: str):
    """
    Determines if a molecule is a cannabinoid based on its SMILES string.
    Cannabinoids typically contain notable features like a resorcinol moiety,
    long hydrocarbon chains often with unsaturations, ether/ester/amide bonds 
    and may feature cyclic structures attached or linked to chains.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as a cannabinoid, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for features that are common in many cannabinoids:
    
    # Check for resorcinol core (common in compounds like THC, CBD)
    resorcinol_pattern = Chem.MolFromSmarts("Oc1ccccc1")
    if mol.HasSubstructMatch(resorcinol_pattern):
        return True, "Classified as cannabinoid by presence of resorcinol core"
    
    # Check for amide & fatty acid connected (typical in certain cannabinoid derivatives)
    amide_fatty_acid_pattern = Chem.MolFromSmarts("C(=O)NCCO")
    if mol.HasSubstructMatch(amide_fatty_acid_pattern):
        return True, "Classified as cannabinoid by presence of an amide bound to ethanolamine"

    # Look for the presence of long hydrocarbon chains with unsaturation (typical in structures like anandamide)
    long_chain_unsaturated_pattern = Chem.MolFromSmarts("CCCC=CCCC=CCCC=CCCC")
    if mol.HasSubstructMatch(long_chain_unsaturated_pattern):
        # Further check for connected ester or ether for specificity
        if any(mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)) for pattern in ["C(=O)OC", "C-O-C"]):
            return True, "Classified as cannabinoid by presence of long unsaturated chain with ester or ether"
        return True, "Classified as cannabinoid by presence of long unsaturated chain"
    
    # Detection for heterocyclic components bracketed to hydrocarbon chains
    cyclic_ether_or_hetero = Chem.MolFromSmarts("c1occc1")
    if mol.HasSubstructMatch(cyclic_ether_or_hetero):
        return True, "Classified as cannabinoid by presence of cyclic ether or other heterocycle"

    # Evaluate using cannabinoid-specific reference framework of known examples
    # If known structural classifications fail, return false
    return False, "Lacks distinct structural features typical of cannabinoids"