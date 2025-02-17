"""
Classifies: CHEBI:23899 icosanoid
"""
"""
Classifies: Icosanoid
Definition: Any member of the group of signalling molecules arising from oxidation of the three C20 essential fatty acids (EFAs) icosapentaenoic acid (EPA), arachidonic acid (AA) and dihomo-gamma-linolenic acid (DGLA).
Heuristic criteria applied:
  • The SMILES must be valid.
  • The molecule should have a total carbon count roughly in the range expected (for the oxidized C20 skeleton, we allow 15–40 carbon atoms).
  • It should contain at least one oxygenated carbonyl function (as in a carboxylic acid or ester group), characteristic in fatty acid derivatives.
  • It should possess at least one or more C=C double bonds (as the EFAs are polyunsaturated).
  • The estimated molecular weight should lie in a typical range (250–900 Da) for icosanoids.
Note: This is a heuristic implementation and may not catch every nuance of icosanoid chemistry.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_icosanoid(smiles: str):
    """
    Determines if a molecule could be an icosanoid based on its SMILES string.
    
    Heuristic criteria:
      - Valid molecule parsing.
      - Carbon count in a reasonable range (e.g., 15 to 40).
      - Contains at least one oxygenated carbonyl (as in acid or ester groups).
      - Contains at least one double bond.
      - Molecular weight falls roughly between 250 and 900 Da.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule passes our heuristic for being an icosanoid, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # Count the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 15 or c_count > 40:
        return False, f"Carbon count {c_count} not in expected range (15-40) for an icosanoid."

    # Compute estimated molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 900:
        return False, f"Molecular weight ({mol_wt:.1f} Da) out of expected range (250-900 Da) for an icosanoid."

    # Check for oxygenated carbonyl group using SMARTS.
    # This pattern looks for a carbonyl (C=O) immediately attached to an oxygen (which could be from a carboxylic acid or ester).
    carbonyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carbonyl_pattern):
        return False, "No carbonyl-oxygen fragment (acid/ester group) found."

    # Count the number of double bonds in the molecule.
    double_bonds = [bond for bond in mol.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE]
    if len(double_bonds) < 1:
        return False, "Too few double bonds; icosanoids are typically derived from polyunsaturated fatty acids."

    # (Optional) Additional check: icosanoids are oxidation products of C20 EFAs.
    # We can check if there is at least one carboxylic acid group.
    acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1,-]")
    if not mol.HasSubstructMatch(acid_pattern):
        # Note: Some icosanoids are esterified or conjugated so they may not show a free acid.
        acid_note = " (No free carboxylic acid detected; molecule may be fully esterified/conjugated.)"
    else:
        acid_note = ""
    
    return True, f"Meets heuristic criteria for an icosanoid. (C_count={c_count}, MW={mol_wt:.1f} Da, double bonds={len(double_bonds)}){acid_note}"

# Example usage:
if __name__ == "__main__":
    # Test a sample icosanoid: leukotriene A4
    smiles_examples = [
        "CCCC\\C=C/C\\C=C/C=C/C=C/[C@@H]1O[C@H]1CCCC(O)=O",  # leukotriene A4
        "CCCC[C@H](OO)\\C=C\\C=C/C\\C=C/C\\C=C/CCCC(O)=O",      # 15(S)-HPETE (modified)
    ]
    for smi in smiles_examples:
        result, reason = is_icosanoid(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")