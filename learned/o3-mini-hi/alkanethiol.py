"""
Classifies: CHEBI:47908 alkanethiol
"""
"""
Classifies: alkanethiol
Definition: An alkanethiol is a compound in which a sulfanyl (-SH) group is attached to an alkyl group.
Examples include 3-Methyl-2-butanethiol, propane-2-thiol, ethanethiol, etc.
This version also excludes molecules that contain amide (peptide) bonds.
"""

from rdkit import Chem

def is_alkanethiol(smiles: str):
    """
    Determines if a molecule is an alkanethiol based on its SMILES string.
    An alkanethiol is defined as a compound with a thiol (–SH) group attached 
    to a saturated alkyl (sp3) carbon. Additionally, molecules that contain amide 
    bonds (e.g. peptide backbones) are rejected.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as an alkanethiol, False otherwise
        str: Explanation for the classification decision
    """
    # Parse the molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Update property cache
    mol.UpdatePropertyCache()
    
    # Exclude molecules that contain amide bonds which are hallmarks of peptides.
    # The SMARTS "C(=O)N" will match substructures having a carbonyl carbon bonded to nitrogen.
    amide_smarts = Chem.MolFromSmarts("C(=O)N")
    if mol.HasSubstructMatch(amide_smarts):
        return False, ("Thiolate group found but molecule contains amide bonds, "
                       "suggesting a peptide or amidic structure rather than a simple alkyl chain.")
    
    valid_thiol_found = False
    invalid_thiol_found = False
    
    # Iterate over atoms looking for sulfur atoms which may represent –SH groups.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 16:  # sulfur
            # Ensure the sulfur has at least one hydrogen (could be implicit)
            if atom.GetTotalNumHs() < 1:
                continue  # not a thiol because no attached hydrogen

            # Look at heavy-atom neighbors only (ignoring hydrogens)
            neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() != 1]
            
            # A genuine thiol typically has S attached to exactly one heavy atom.
            if len(neighbors) != 1:
                invalid_thiol_found = True
                continue

            alkyl_candidate = neighbors[0]
            
            # The neighbor should be a carbon
            if alkyl_candidate.GetAtomicNum() != 6:
                invalid_thiol_found = True
                continue

            # The attached carbon must be saturated (sp3 hybridized) and not aromatic.
            if (alkyl_candidate.GetHybridization() != Chem.rdchem.HybridizationType.SP3 or 
                alkyl_candidate.GetIsAromatic()):
                invalid_thiol_found = True
                continue

            # We have found a valid –SH attached to a saturated alkyl group.
            valid_thiol_found = True
            break

    if valid_thiol_found:
        return True, "Contains at least one –SH group attached to a saturated alkyl (sp3 carbon) chain."
    elif invalid_thiol_found:
        return False, "Thiolate group found but it is not attached to a saturated alkyl (sp3 carbon) chain."
    else:
        return False, "No thiol (–SH) group found in the molecule."

# Example usage when running the script directly:
if __name__ == "__main__":
    test_smiles = [
        "SC(C(C)C)C",         # 3-Methyl-2-butanethiol (true positive)
        "CC(C)S",             # propane-2-thiol (true positive)
        "SC(C(C)C)CO",        # xi-2-Mercapto-3-methyl-1-butanol (true positive)
        "S\C=C(\CC)/C",       # 2-Methyl-1-propenethiol (false; unsaturated carbon)
        "SC[C@H](NC(=O)[C@@H](N)CCCN=C(N)N)C(=O)N"  # part of a peptide (false positive originally)
    ]
    
    for smi in test_smiles:
        is_valid, reason = is_alkanethiol(smi)
        print(f"SMILES: {smi:50s} -> {is_valid} ({reason})")