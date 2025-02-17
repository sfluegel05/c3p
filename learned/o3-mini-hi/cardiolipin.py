"""
Classifies: CHEBI:28494 cardiolipin
"""
#!/usr/bin/env python
"""
Classifies: cardiolipin
A phosphatidylglycerol composed of two molecules of phosphatidic acid covalently linked to a molecule of glycerol.
Heuristics:
  - Must have exactly 2 phosphorus atoms.
  - Must contain at least 4 acyl ester bonds (OC(=O)) for 4 fatty acyl chains.
  - Must include a glycerol linker in the headgroup.
"""

from rdkit import Chem

def is_cardiolipin(smiles: str):
    """
    Determines if a molecule is a cardiolipin based on its SMILES string.
    Cardiolipin is defined as a phosphatidylglycerol composed of two molecules of phosphatidic acid
    covalently linked via a glycerol linker. In other words, cardiolipin should have two phosphate groups,
    four acyl ester bonds typical of fatty acyl chains, and a glycerol linking motif in its headgroup.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a cardiolipin, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for phosphorus atoms - cardiolipin should have exactly 2 phosphate groups.
    phosphorus_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if len(phosphorus_atoms) != 2:
        return False, f"Expected 2 phosphorus atoms, found {len(phosphorus_atoms)}"
    
    # Look for acyl ester bonds "OC(=O)" that indicate fatty acyl chains.
    ester_pattern = Chem.MolFromSmarts("OC(=O)")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 4:
        return False, f"Expected at least 4 acyl ester groups, found {len(ester_matches)}"

    # Look for a glycerol linker in the headgroup.
    # Instead of our previous strict pattern, we now use a more general pattern.
    # The idea here is to capture a fragment like "OCC(O)OP" that is present in many valid cardiolipin headgroups.
    glycerol_linker_pattern = Chem.MolFromSmarts("OCC(O)OP")
    if not mol.HasSubstructMatch(glycerol_linker_pattern):
        return False, "Glycerol linker not detected in the headgroup"
    
    # If all heuristics pass, classify as cardiolipin.
    return True, "Contains 2 phosphorus atoms, at least 4 acyl ester groups, and a glycerol linker consistent with cardiolipin"

# Example usage (for testing purposes)
if __name__ == "__main__":
    # Example SMILES for cardiolipin from provided examples.
    example_smiles = "P(OC[C@H](OC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCC/C=C\\CCCCCCCC)(OC[C@H](O)COP(OC[C@H](OC(=O)CCCCCCC/C=C\\C/C=C\\CCCCC)COC(=O)CCCCCCC/C=C\\C/C=C\\CCCCC)(O)=O)(O)=O"
    result, reason = is_cardiolipin(example_smiles)
    print("Is cardiolipin?", result)
    print("Reason:", reason)