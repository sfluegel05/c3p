"""
Classifies: CHEBI:28765 phosphatidylinositol phosphate
"""
"""
Classifies: Phosphatidylinositol phosphate (a phosphoinositide, one of seven naturally occurring).
A phosphatidylinositol phosphate is expected to contain:
  - A (myo-)inositol head group (a six-membered ring with several hydroxyl groups),
  - A phosphate ester linking the inositol to a glycerol backbone,
  - At least two long-chain acyl groups attached as esters (forming the diacylglycerol part).

This program is an approximate classifier that uses chemical substructure matching via RDKit.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylinositol_phosphate(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol phosphate based on its SMILES string.
    The classification is based on the presence of a (myo-)inositol ring, a phosphate linking group,
    and at least two acyl chain ester groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a phosphatidylinositol phosphate, False otherwise.
        str: Reason for classification
    """
    # Parse the SMILES to an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check for an inositol head group.
    # This SMARTS looks for a six-membered carbon ring with multiple hydroxyl (O) substituents.
    # (This is an approximate pattern for myo-inositol.)
    inositol_smarts = "C1(C(C(C(C(C1O)O)O)O)O)"  # non-stereo version – will match several inositol forms.
    inositol_pattern = Chem.MolFromSmarts(inositol_smarts)
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "Inositol head group not found"
    
    # 2. Check for phosphate ester functionality.
    # Look for a phosphate group linked via an oxygen.
    phosphate_smarts = "OP(=O)(O)O"
    phosphate_pattern = Chem.MolFromSmarts(phosphate_smarts)
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "Phosphate ester group not found"
    
    # 3. Check for diacylglycerol: typically at least two acyl chains attached via ester bonds.
    # Here we search for ester fragments that resemble a fatty acid chain.
    # This SMARTS looks for an ester: an oxygen attached to a C=O which is in turn connected to a chain of at least 6 carbons.
    acyl_smarts = "OC(=O)CCCCCC"  
    acyl_pattern = Chem.MolFromSmarts(acyl_smarts)
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if len(acyl_matches) < 2:
        return False, f"Found only {len(acyl_matches)} acyl ester group(s), need at least 2"
    
    # 4. As an extra check: look for the presence of phosphorus (should be at least one).
    p_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if len(p_atoms) < 1:
        return False, "No phosphorus atom found"
    
    # Optionally, one may wish to check overall molecular size. Phosphatidylinositol phosphates are generally large.
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 500:
        return False, f"Molecular weight ({mw:.1f} Da) is lower than expected for a phosphatidylinositol phosphate"
    
    # If all conditions pass, then we classify the compound as a phosphatidylinositol phosphate.
    return True, "Molecule contains the inositol head group, a linking phosphate ester, and at least two acyl chains"

# Example usage:
if __name__ == "__main__":
    # Using one of the provided example SMILES:
    test_smiles = "CCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@H]1[C@H](O)[C@@H](OP(O)(O)=O)[C@H](OP(O)(O)=O)[C@H]1O)OC(=O)CCCCCCC"
    result, reason = is_phosphatidylinositol_phosphate(test_smiles)
    print(result, reason)