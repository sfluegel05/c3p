"""
Classifies: CHEBI:67194 cannabinoid
"""
"""
Classifies: cannabinoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_cannabinoid(smiles: str):
    """
    Determines if a molecule is a cannabinoid based on its SMILES string.
    Cannabinoids are classified based on their core structures such as:
    - Phytocannabinoids like THC, CBD with a dibenzopyran ring system.
    - Endocannabinoids like anandamide, which are N-acylethanolamines of long-chain fatty acids.
    - Synthetic cannabinoids with specific core structures like indoles or indazoles with alkyl side chains.

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

    # Define patterns

    # Pattern for THC-like cannabinoids (dibenzopyran ring system)
    thc_pattern = Chem.MolFromSmarts("c1cc2c(c1)c(O)cc(C(C)(C)O)c2O")
    if mol.HasSubstructMatch(thc_pattern):
        return True, "Contains dibenzopyran core like THC"

    # Pattern for CBD-like cannabinoids (resorcinol moiety with alkyl chain)
    cbd_pattern = Chem.MolFromSmarts("c1cc(cc(c1O)O)CC=C(C)CCC=C(C)C")
    if mol.HasSubstructMatch(cbd_pattern):
        return True, "Contains resorcinol moiety with alkyl chain like CBD"

    # Pattern for anandamide-like endocannabinoids (N-acylethanolamine of arachidonic acid)
    anandamide_pattern = Chem.MolFromSmarts("CCCCC=CCCC=CCCC=CCCC(=O)NCCO")
    if mol.HasSubstructMatch(anandamide_pattern):
        return True, "Contains N-acylethanolamine of long-chain fatty acid like anandamide"

    # Pattern for 2-arachidonoylglycerol-like molecules (glycerol ester of arachidonic acid)
    glycerol_ester_pattern = Chem.MolFromSmarts("OCC(O)COC(=O)CCCCC=CCCC=CCCC=CCC")
    if mol.HasSubstructMatch(glycerol_ester_pattern):
        return True, "Contains glycerol ester of arachidonic acid"

    # Pattern for synthetic cannabinoids with indole core and alkyl side chain
    indole_pattern = Chem.MolFromSmarts("c1ccc2c(c1)[nH]cc2C(=O)C3=CC=CC(I)=C3")
    if mol.HasSubstructMatch(indole_pattern):
        return True, "Contains indole core with specific substitution like synthetic cannabinoids"

    # Pattern for synthetic cannabinoids with indazole core
    indazole_pattern = Chem.MolFromSmarts("c1ccc2c(c1)nnc2C(=O)N3CCCCC3")
    if mol.HasSubstructMatch(indazole_pattern):
        return True, "Contains indazole core with specific substitution like synthetic cannabinoids"

    # Pattern for oxygen-containing heterocyclic ring connected to alkyl chain
    heterocycle_alkyl_pattern = Chem.MolFromSmarts("c1oc(cc1)C(C)=C")
    if mol.HasSubstructMatch(heterocycle_alkyl_pattern):
        # Check for long alkyl chain
        alkyl_chain_pattern = Chem.MolFromSmarts("CCCCCCCC")
        if mol.HasSubstructMatch(alkyl_chain_pattern):
            return True, "Contains oxygen heterocycle connected to long alkyl chain"

    # Pattern for N-acylethanolamines with long-chain fatty acids
    n_acylethanolamine_pattern = Chem.MolFromSmarts("C(=O)NCCO")
    if mol.HasSubstructMatch(n_acylethanolamine_pattern):
        # Check for long-chain fatty acid
        fatty_acid_chain = Chem.MolFromSmarts("CCCCCCCCCCCCCCCCCC")
        if mol.HasSubstructMatch(fatty_acid_chain):
            return True, "Contains N-acylethanolamine of long-chain fatty acid"

    # If none of the specific patterns match, it's not classified as a cannabinoid
    return False, "Does not match characteristic structures of cannabinoids"