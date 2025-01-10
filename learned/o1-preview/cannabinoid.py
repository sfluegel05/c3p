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
    Cannabinoids are classified into:
    - Phytocannabinoids: compounds like THC and CBD, containing a dibenzopyran or similar core structure.
    - Endocannabinoids: endogenous lipid-based neurotransmitters like anandamide and 2-AG, derivatives of arachidonic acid.
    - Synthetic cannabinoids: compounds often containing indole or indazole core structures with specific substitutions.

    This function checks for these characteristic structures.

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

    # Pattern for phytocannabinoids (e.g., THC, CBD) - dibenzopyran or resorcinol core
    phytocannabinoid_pattern = Chem.MolFromSmarts("c1cc2c(cc1)c(O)cc(O)c2")
    # This matches the resorcinol moiety fused with a benzene ring

    if mol.HasSubstructMatch(phytocannabinoid_pattern):
        return True, "Contains dibenzopyran or resorcinol core like phytocannabinoids (e.g., THC, CBD)"

    # Pattern for anandamide (N-arachidonoylethanolamine)
    # N-acylethanolamine with 20-carbon chain and 4 cis double bonds
    anandamide_pattern = Chem.MolFromSmarts("OCCNC(=O)CCCCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CC")
    if mol.HasSubstructMatch(anandamide_pattern):
        return True, "Contains N-arachidonoylethanolamine structure like anandamide"

    # Pattern for 2-arachidonoylglycerol (2-AG)
    # Monoacylglycerol with arachidonoyl chain
    monoacylglycerol_pattern = Chem.MolFromSmarts("OCC(O)COC(=O)CCCCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CC")
    if mol.HasSubstructMatch(monoacylglycerol_pattern):
        return True, "Contains 2-arachidonoylglycerol structure"

    # Pattern for synthetic cannabinoids with indole core and alkyl side chain
    indole_pattern = Chem.MolFromSmarts("c1ccc2c(c1)[nH]cc2")
    alkyl_chain_pattern = Chem.MolFromSmarts("C[CH2][CH2][CH2][CH2][CH2]")  # pentyl chain
    if mol.HasSubstructMatch(indole_pattern) and mol.HasSubstructMatch(alkyl_chain_pattern):
        return True, "Contains indole core with alkyl side chain like synthetic cannabinoids"

    # Pattern for synthetic cannabinoids with indazole core
    indazole_pattern = Chem.MolFromSmarts("c1ccc2[nH]ncc2c1")
    if mol.HasSubstructMatch(indazole_pattern) and mol.HasSubstructMatch(alkyl_chain_pattern):
        return True, "Contains indazole core with alkyl side chain like synthetic cannabinoids"

    # Pattern for cannabinoids derived from arachidonic acid
    arachidonic_acid_pattern = Chem.MolFromSmarts("CCCCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CC")
    arachidonic_derivative_pattern = Chem.MolFromSmarts("C(=O)O" + arachidonic_acid_pattern)  # Ester or acid derivatives
    if mol.HasSubstructMatch(arachidonic_derivative_pattern):
        return True, "Contains arachidonic acid derivative"

    # If none of the specific patterns match, it's not classified as a cannabinoid
    return False, "Does not match characteristic structures of cannabinoids"