"""
Classifies: CHEBI:28494 cardiolipin
"""
#!/usr/bin/env python
"""
Classifies: cardiolipin
A phosphatidylglycerol composed of two molecules of phosphatidic acid covalently linked to a molecule of glycerol.
Heuristics:
  - Must have exactly 2 phosphorus atoms.
  - Must contain 4 acyl ester bonds (OC(=O)) corresponding to 4 fatty acyl chains.
  - Must include at least one glycerol unit in the headgroup.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cardiolipin(smiles: str):
    """
    Determines if a molecule is a cardiolipin based on its SMILES string.
    Cardiolipin is defined as a phosphatidylglycerol composed of two molecules of phosphatidic acid
    covalently linked to a molecule of glycerol. In other words, cardiolipin has two phosphate groups,
    four acyl ester bonds and a glycerol linker in its headgroup.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as a cardiolipin, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for phosphorus atoms – cardiolipin should have exactly 2 phosphate groups.
    phosphorus_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if len(phosphorus_atoms) != 2:
        return False, f"Expected 2 phosphorus atoms, found {len(phosphorus_atoms)}"

    # Look for acyl ester bonds typical for fatty acid chains.
    # Here our ester SMARTS "OC(=O)" should capture an oxygen linked to a carbonyl,
    # representing the ester bond from a glycerol to a fatty acyl chain.
    ester_pattern = Chem.MolFromSmarts("OC(=O)")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 4:
        return False, f"Expected at least 4 acyl ester groups, found {len(ester_matches)}"

    # Look for a glycerol unit as a bridging (linker) moiety.
    # This SMARTS is a heuristic for a glycerol backbone fragment: CH2-O-CH(OH)-CH2.
    glycerol_pattern = Chem.MolFromSmarts("[CH2][O][CH]([OH])[CH2]")
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    if len(glycerol_matches) < 1:
        return False, "Glycerol linker not detected in the headgroup"

    # Optionally, one could include further checks on the number of fatty acyl chains (by, e.g., evaluating the length
    # of hydrocarbon chains, rotatable bonds, or molecular weight) but here we use only substructure matches.

    return True, "Contains 2 phosphorus atoms, 4 acyl ester groups, and a glycerol linker consistent with cardiolipin"

# Example usage (for testing purposes; remove or modify when integrating as a module)
if __name__ == "__main__":
    # One example SMILES for cardiolipin (from the provided examples)
    example_smiles = "P(OC[C@H](OC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCC/C=C\\CCCCCCCC)(OC[C@H](O)COP(OC[C@H](OC(=O)CCCCCCC/C=C\\C/C=C\\CCCCC)COC(=O)CCCCCCC/C=C\\C/C=C\\CCCCC)(O)=O)(O)=O"
    result, reason = is_cardiolipin(example_smiles)
    print("Is cardiolipin?", result)
    print("Reason:", reason)